/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mpipar;

import cetus.analysis.DDTDriver;
import cetus.hir.Annotation;
import cetus.hir.ArraySpecifier;
import cetus.hir.AssignmentExpression;
import cetus.hir.AssignmentOperator;
import cetus.hir.BinaryExpression;
import cetus.hir.BinaryOperator;
import cetus.hir.CetusAnnotation;
import cetus.hir.CommentAnnotation;
import cetus.hir.CompoundStatement;
import cetus.hir.Declarator;
import cetus.hir.DepthFirstIterator;
import cetus.hir.Expression;
import cetus.hir.ExpressionStatement;
import cetus.hir.FlatIterator;
import cetus.hir.ForLoop;
import cetus.hir.FunctionCall;
import cetus.hir.IDExpression;
import cetus.hir.IfStatement;
import cetus.hir.IntegerLiteral;
import cetus.hir.Loop;
import cetus.hir.NameID;
import cetus.hir.OmpAnnotation;
import cetus.hir.PragmaAnnotation;
import cetus.hir.Procedure;
import cetus.hir.Program;
import cetus.hir.Specifier;
import cetus.hir.Statement;
import cetus.hir.SymbolTools;
import cetus.hir.Symbolic;
import cetus.hir.TranslationUnit;
import cetus.hir.Traversable;
import cetus.hir.UnaryExpression;
import cetus.hir.UnaryOperator;
import cetus.hir.VariableDeclaration;
import cetus.hir.VariableDeclarator;
import cetus.transforms.LoopTransformPass;
import java.util.ArrayList;
import java.util.List;

/**
 * parse MPI parallel - MPI accessedData - MPI reduction - annotations and find
 * the neccessary infos to partition parallel loops - and add communication
 * codes
 *
 * @author tibebu
 */
public class MPI_LoopTransform extends LoopTransformPass {

    public String comm_scheme_selected;
     /**
     * with the broadcast communication from master to slaves - in the 2nd comm
     * the default 
     */
    public static String COMM_SCHEME_ONE = "comm_scheme_one";
    /**
     * one to one communication - from master to slave - in 2nd comm
     * communicates elements of the data except those modified by the slave 
     */
    public static String COMM_SCHEME_TWO = "comm_scheme_two";
    
    public boolean timeInstrumenting;
    public MPI_LoopTransform(Program program,String commScheme,boolean timeInstrumenting) {
        super(program);
        this.comm_scheme_selected=commScheme;
        this.timeInstrumenting = timeInstrumenting;
        //disable analysis consistency checking
        disable_protection = true;
        disable_invalidation = true;
        //if time instrumenting
        if(this.timeInstrumenting)
            insertTimeInstrumentation(program);
    }

    @Override
    public void transformLoop(Loop arg0) {
        if (arg0 instanceof ForLoop) {
            ForLoop floop = (ForLoop) arg0;
          
            //remove notations other than the MPI_Annotation
            deatachNotationsExceptMPIAnnotations(floop);
            //Now try Mpi parallelizing MPI parallel loops - first check if it is mpi parallel & accessedvariables are found
            if (floop.containsAnnotation(MPI_Annotation.class, MPI_Annotation.KEY_PARALLEL) && floop.containsAnnotation(MPI_Annotation.class, MPI_Annotation.KEY_ACCESSED_DATA)) {
                //get the depth - 
                int depth = ((MPI_Annotation) floop.getAnnotation(MPI_Annotation.class, MPI_Annotation.KEY_PARALLEL)).get(MPI_Annotation.KEY_PARALLEL);
                //if depth - outer paralle loop 
                //only parallelize outer parallel loops
                //parse and get shared variables
                if (depth == 0) {
                    MPI_Annotation mpiSVannot = floop.getAnnotation(MPI_Annotation.class, MPI_Annotation.KEY_ACCESSED_DATA);
                    List<MPI_AccessedData> accessedVariables = mpiSVannot.get(MPI_Annotation.KEY_ACCESSED_DATA);

                    if (!accessedVariables.isEmpty()) {
                        //insert computation partitioning and mapping code
                        //loopComputationPartition(floop);
                        insertComputationPartitioningMappingCode(floop);
                        //communication code generation - use default comm scheme
                        insertCommunicationCodes(floop, accessedVariables,this.comm_scheme_selected,this.timeInstrumenting);

                    }
                }

            }

            //after transformation  regenerate loop dependence
            //DDTDriver ddt = new DDTDriver(program);
            // ddt.start();
        }
    }

    @Override
    public String getPassName() {
        return "MPI Loop trasnform - computation partitioning";
    }

    /**
     * remove annotation except the MPI_Annotations
     *
     * @param loop
     */
    public void deatachNotationsExceptMPIAnnotations(ForLoop loop) {
        //get the MPIAnnotations first
        MPI_Annotation mpiAnnotationPar = loop.getAnnotation(MPI_Annotation.class, MPI_Annotation.KEY_PARALLEL);
        MPI_Annotation mpiAnnotationSv = loop.getAnnotation(MPI_Annotation.class, MPI_Annotation.KEY_ACCESSED_DATA);
        MPI_Annotation mpiAnnotationRv = loop.getAnnotation(MPI_Annotation.class, MPI_Annotation.KEY_REDUCTION);
        MPI_Annotation mpiAnnotationNsc = loop.getAnnotation(MPI_Annotation.class, MPI_Annotation.KEY_NO_SECOND_COMM);
        

        loop.removeAnnotations();
        if (mpiAnnotationPar != null) {
            loop.annotate(mpiAnnotationPar);
        }
        if (mpiAnnotationRv != null) {
            loop.annotate(mpiAnnotationRv);
        }
        if (mpiAnnotationSv != null) {
            loop.annotate(mpiAnnotationSv);
        }
        if (mpiAnnotationNsc != null) {
            loop.annotate(mpiAnnotationNsc);
        }

    }

   
    /**
     * inserts loop partitioning and mapping code
     * @param loop 
     */
    public void insertComputationPartitioningMappingCode(ForLoop loop){
        CompoundStatement container = new CompoundStatement(); //container for partitioning code
        container.annotate(new CommentAnnotation("Computation partitioning and mapping code begins "));
        //1- declare partitioning variables (int Ic,Rc,Lc,Uc,Lcp,Ucp)- scope loop
        VariableDeclaration Ic = MPI_Utility.DeclareVariable(container, "Ic", Specifier.INT);
        VariableDeclaration Rc = MPI_Utility.DeclareVariable(container, "Rc", Specifier.INT);
        VariableDeclaration Lc = MPI_Utility.DeclareVariable(container, "Lc", Specifier.INT);
        VariableDeclaration Uc = MPI_Utility.DeclareVariable(container, "Uc", Specifier.INT);
        VariableDeclaration Lcp = MPI_Utility.DeclareVariable(container, "Lcp", Specifier.INT);
        VariableDeclaration Ucp = MPI_Utility.DeclareVariable(container, "Ucp", Specifier.INT);
        //2- get loop bounds [Lc,Uc) // 
        Statement stmntLc = new ExpressionStatement(new AssignmentExpression(Lc.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, MPI_Utility.findLoopLowerBound(loop).clone()));
        Statement stmntUc = new ExpressionStatement(new AssignmentExpression(Uc.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, MPI_Utility.FindLoopUpperBound(loop).clone()));
        container.addStatement(stmntLc);
        container.addStatement(stmntUc);
        //3- compute size of computation and computation partition remains
        //Ic = Uc-Lc; Rc = Ic%size;
        Statement stmntIc = new ExpressionStatement(new AssignmentExpression(Ic.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, new BinaryExpression(Uc.getDeclaredIDs().get(0).clone(), BinaryOperator.SUBTRACT, Lc.getDeclaredIDs().get(0).clone())));
        Statement stmntRc = new ExpressionStatement(new AssignmentExpression(Rc.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, new BinaryExpression(Ic.getDeclaredIDs().get(0).clone(), BinaryOperator.MODULUS, new NameID("size"))));
        container.addStatement(stmntIc);
        container.addStatement(stmntRc);
        //4 compute computation partition bounds [Lcp,Ucp)
        //if(rank==0){Lcp=Lc;Ucp=Lcp+Ic/size+Rc}else{Lcp=Lc+rank*Ic/size+Rc);Ucp=Lcp+Ic/size}
        IfStatement ccpb = new IfStatement(new BinaryExpression(new NameID("rank"), BinaryOperator.COMPARE_EQ, new IntegerLiteral(0)), new CompoundStatement());
        CompoundStatement ccpb_true_clause = new CompoundStatement();
        CompoundStatement ccpb_false_clause = new CompoundStatement();
        ccpb.setThenStatement(ccpb_true_clause);
        ccpb.setElseStatement(ccpb_false_clause);
        //for the true clause - master node
        Statement stmntLcp_0 = new ExpressionStatement(new AssignmentExpression(Lcp.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, Lc.getDeclaredIDs().get(0).clone()));
            //Ic_size = 'Ic/size'
        Expression Ic_size = new BinaryExpression(Ic.getDeclaredIDs().get(0).clone(), BinaryOperator.DIVIDE, new NameID("size"));
        Statement stmntUcp_0 = new ExpressionStatement(new AssignmentExpression(Ucp.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, new BinaryExpression(new BinaryExpression(Lcp.getDeclaredIDs().get(0).clone(),BinaryOperator.ADD,Ic_size.clone()), BinaryOperator.ADD,Rc.getDeclaredIDs().get(0).clone() )));
        ccpb_true_clause.addStatement(stmntLcp_0);
        ccpb_true_clause.addStatement(stmntUcp_0);
        //for the false clause - other than master node
        //rank_Ic_size = rank*Ic/size
        Expression rank_Ic_size = new BinaryExpression(new NameID("rank"), BinaryOperator.MULTIPLY, Ic_size.clone());
        Statement stmntLcp_1 = new ExpressionStatement(new AssignmentExpression(Lcp.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL,new BinaryExpression(new BinaryExpression(Lc.getDeclaredIDs().get(0).clone(),BinaryOperator.ADD,rank_Ic_size.clone()), BinaryOperator.ADD,Rc.getDeclaredIDs().get(0).clone())));
        Statement stmntUcp_1 = new ExpressionStatement(new AssignmentExpression(Ucp.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL,new BinaryExpression(Lcp.getDeclaredIDs().get(0).clone(), BinaryOperator.ADD, Ic_size.clone())));
        ccpb_false_clause.addStatement(stmntLcp_1);
        ccpb_false_clause.addStatement(stmntUcp_1);
        container.addStatement(ccpb);
        //5 computation partition mapping  cpm
        //if(index>=Lcp && index<Ucp){loop.body}
        Expression loopIndex = MPI_Utility.FindLoopIndex(loop).clone();
        //index>=Lcp
        Expression index_lcp = new BinaryExpression(loopIndex.clone(), BinaryOperator.COMPARE_GE, Lcp.getDeclaredIDs().get(0).clone());
        //index<Ucp
        Expression index_Ucp = new BinaryExpression(loopIndex.clone(), BinaryOperator.COMPARE_LT, Ucp.getDeclaredIDs().get(0).clone());
        IfStatement cpm = new IfStatement(new BinaryExpression(index_lcp, BinaryOperator.LOGICAL_AND, index_Ucp), new CompoundStatement());
        cpm.setThenStatement(loop.getBody().clone());
        container.addStatement(cpm);
        
        //set loop body
        loop.setBody(container);
        
    }
    /**
     * insert communication codes 
     * @param floop
     * @param accessedVariables 
     */
    public void insertCommunicationCodes(ForLoop floop, List<MPI_AccessedData> accessedVariables,String comm_scheme, boolean timeReport) {
        List<MPI_AccessedData> notReductionsharedVariables = new ArrayList<>();
        List<MPI_AccessedData> reductionsharedVariables = new ArrayList<>();
        for (MPI_AccessedData sharedVariable : accessedVariables) {
            if (sharedVariable.isIsReductionVariable()) {
                reductionsharedVariables.add(sharedVariable);
            } else {
                notReductionsharedVariables.add(sharedVariable);
            }
        }

        //generate the communication code for shared variables that are not reduction variables & 
        if (!notReductionsharedVariables.isEmpty()) {
            //CommunicationCodeGeneration4NonReductionAccesssedData(floop, notReductionsharedVariables,comm_scheme);
            communicationCodeGeneration4NonReductionAccessedData(floop, notReductionsharedVariables,timeReport);
        }
        //generate comm code gen for shared variables - thatare reductionVariables and lastindexSubscript is not function of the paralle loop index    
        if (!reductionsharedVariables.isEmpty()) {
            CommunicationCodeGeneration4ReductionSharedVariables(floop, reductionsharedVariables);
        }
    }

    /**
     * inserts communication code to accesseddata that are reduction variables
     * @param loop
     * @param sharedVariables 
     */
    public void CommunicationCodeGeneration4ReductionSharedVariables(ForLoop loop, List<MPI_AccessedData> sharedVariables) {
        CompoundStatement CommunicationCodeGeneration4ReductionSharedVariables = new CompoundStatement();
        CommunicationCodeGeneration4ReductionSharedVariables.annotate(new CommentAnnotation("COMMUNICATION CODE FOR SHARED - REDUCTION VARIABLES BEGIN"));
        //for each shared variables
        for (MPI_AccessedData sharedVariable : sharedVariables) {

            //declare the variable to store the result- variable has datatype of the sv, & its scope will be the loop parent
            VariableDeclaration tempResult = MPI_Utility.DeclareVariable(CommunicationCodeGeneration4ReductionSharedVariables, "tempResult", sharedVariable.getDataType());

            tempResult.annotate(new CommentAnnotation("For Reduction variable -" + sharedVariable.getVariableSymbol().getSymbolName()));
            //add the reduce mpi function
            //MPI_Reduce(&C[i][j],&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            FunctionCall MPI_Reduce = new FunctionCall(SymbolTools.getOrphanID("MPI_Reduce"));
            //add the argument - in variable
            String argument1 = "&" + sharedVariable.getVarExpression(); //if the variable is array
            if (!sharedVariable.isArray()) //if it is scalar
            {
                argument1 = "&" + sharedVariable.getVariableSymbol().getSymbolName();
            }
            MPI_Reduce.addArgument(new NameID(argument1));
            //add the out variable argument - temp buffer to store the collected data
            String argument2 = "&" + tempResult.getDeclaredIDs().get(0).clone();
            MPI_Reduce.addArgument(new NameID(argument2));
            //add the 3rd argument - count - num of data
            MPI_Reduce.addArgument(new IntegerLiteral(1));
            //add the 4th argument - datatype
            String dataType = "MPI_" + sharedVariable.getDataType().toString().toUpperCase();
            MPI_Reduce.addArgument(new NameID(dataType));
            //add the 5th argum - the association/ reduction operation
            String reductionOp = "";
            if (sharedVariable.getReductionOperator().equals(BinaryOperator.ADD)) {
                reductionOp = "MPI_SUM";
            } else if (sharedVariable.getReductionOperator().equals(BinaryOperator.MULTIPLY)) {
                reductionOp = "MPI_PROD";
            }
            //TODO - include the rest reduction operations
            //MPI_MAX, MPI_MIN,          
            MPI_Reduce.addArgument(new NameID(reductionOp));
            //add the 6th argument - destination - master node
            MPI_Reduce.addArgument(new IntegerLiteral(0));
            //add the final argument
            MPI_Reduce.addArgument(new NameID("MPI_COMM_WORLD"));
            Statement MPI_ReduceStat = new ExpressionStatement(MPI_Reduce);
            //put the result to the shared variable - only the master do this
            //if(rank==0){ C[i][j] = result}
            Expression ifStatCondition = new BinaryExpression(new NameID("rank"), BinaryOperator.COMPARE_EQ, new IntegerLiteral(0));
            CompoundStatement trueCalsue = new CompoundStatement();
            IfStatement ifsSta = new IfStatement(ifStatCondition, trueCalsue);
            //sv=result
            Expression varExp = null;
            if (sharedVariable.isArray()) {
                varExp = sharedVariable.getVarExpression().clone();
            } else {
                varExp = new NameID(sharedVariable.getVariableSymbol().getSymbolName());
            }
            Statement st = new ExpressionStatement(new AssignmentExpression(varExp.clone(), AssignmentOperator.NORMAL, tempResult.getDeclaredIDs().get(0).clone()));
            trueCalsue.addStatement(st);

            //add the reduce comm to compound statment
            CommunicationCodeGeneration4ReductionSharedVariables.addStatement(MPI_ReduceStat);
            CommunicationCodeGeneration4ReductionSharedVariables.addStatement(ifsSta);

        }

        //now add the statments to the loop parent -- after the loop
        CompoundStatement par = (CompoundStatement) loop.getParent();
        par.addStatementAfter(loop, CommunicationCodeGeneration4ReductionSharedVariables);

    }

    /**
     * insert brodcasting and/or gathering comm code - after the loop - for each accessedData
     * @param loop
     * @param accessedDataList 
     */
    public void communicationCodeGeneration4NonReductionAccessedData(ForLoop loop,List<MPI_AccessedData> accessedDataList,boolean timeReport){
        //step 1 - create the 'if(size>1){}' container
        CompoundStatement container = new CompoundStatement();
        
        //insert time instrumenting - start mpipar_comm_timer
        if(timeReport){
             FunctionCall start_timer=  new FunctionCall(SymbolTools.getOrphanID("mpipar_comm_timer_start"));
                  //prepare mpipar_comm_timer_start(); statement
            Statement time_start_stmnt = new ExpressionStatement(start_timer);
            container.addStatement(time_start_stmnt);
        }
      
        
        IfStatement ifContainer = new IfStatement(new BinaryExpression(new NameID("size"), BinaryOperator.COMPARE_GT, new IntegerLiteral(1)), container);
        ifContainer.annotate(new CommentAnnotation("COMMUNICATION CODE STARTS"));
        //step 2 - declare comm variables Ud,Ld,Id,Rd,CD_L_p,CD_U_p,CD_S_p,PSM
        VariableDeclaration Ld = MPI_Utility.DeclareVariable(container, "Ld", Specifier.INT);
        VariableDeclaration Ud = MPI_Utility.DeclareVariable(container, "Ud", Specifier.INT);
        VariableDeclaration Id = MPI_Utility.DeclareVariable(container, "Id", Specifier.INT);
        VariableDeclaration Rd = MPI_Utility.DeclareVariable(container, "Rd", Specifier.INT);
        VariableDeclaration CD_L_p = MPI_Utility.DeclareVariable(container, "CD_L_p", Specifier.INT);
        VariableDeclaration CD_U_p = MPI_Utility.DeclareVariable(container, "CD_U_p", Specifier.INT);
        VariableDeclaration CD_S_p = MPI_Utility.DeclareVariable(container, "CD_S_p", Specifier.INT);
        VariableDeclaration PSM = MPI_Utility.DeclareVariable(container, "PSM", Specifier.INT);
        
        //get gathering and broadcasting accessedData List
        List<MPI_AccessedData> gatherADL = new ArrayList<>();
        List<MPI_AccessedData> broadcastADL = new ArrayList<>();
        for (MPI_AccessedData accessedData1 : accessedDataList) {
            if(accessedData1.getCommType().equals(MPI_AccessedData.COMMUNICATION_TYPE_BROADCASTING))
                broadcastADL.add(accessedData1);
            else
                gatherADL.add(accessedData1);
        }
        //step 3 comm code for broadcastADL
        if(!broadcastADL.isEmpty()){
            //declare forloop index , int p
            VariableDeclaration p = MPI_Utility.DeclareVariable(container, "p", Specifier.INT);
         
            ForLoop broadcastLoop = new ForLoop(null, null, null, null);
             //set the loop initial statment // p=0 ; set the condition // p<size  // set the step p++
            broadcastLoop.setInitialStatement(new ExpressionStatement(new AssignmentExpression(p.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, new IntegerLiteral(0))));
            broadcastLoop.setCondition(new BinaryExpression(p.getDeclaredIDs().get(0).clone(), BinaryOperator.COMPARE_LT, new NameID("size")));
            broadcastLoop.setStep(new UnaryExpression(UnaryOperator.POST_INCREMENT, p.getDeclaredIDs().get(0).clone()));
            //declare the forloop body 
            CompoundStatement broadcastLoopBody = new CompoundStatement();
            broadcastLoop.setBody(broadcastLoopBody);
            //iterate throgh each broadcastADL
            for (MPI_AccessedData broadcastADL1 : broadcastADL) {
                //creat container
                CompoundStatement badContainer = new CompoundStatement();
                badContainer.annotate(new CommentAnnotation("For accessed data - "+broadcastADL1.getVariableSymbol().getSymbolName()+" - all-to-all broadcasting"));
                //step 3.1 add statments - that compute comm variables
                //first find the CD partitioned dimension rangeAccess/bounds lb & ub
                Expression lb = broadcastADL1.getAccessRange(broadcastADL1.getDimOfIndexOftheParallelLoop()).getLB().clone();
                Expression ub = broadcastADL1.getAccessRange(broadcastADL1.getDimOfIndexOftheParallelLoop()).getUB().clone();
                
                Statement stmntLd = new ExpressionStatement(new AssignmentExpression(Ld.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, lb.clone()));
                Statement stmntUd = new ExpressionStatement(new AssignmentExpression(Ud.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, ub.clone()));
                badContainer.addStatement(stmntLd);
                badContainer.addStatement(stmntUd);
                //a- compute size of coomunication data and data partition remains
                //Id = Ud-Ld; Rd = Id%size;
                Statement stmntId = new ExpressionStatement(new AssignmentExpression(Id.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, new BinaryExpression(Ud.getDeclaredIDs().get(0).clone(), BinaryOperator.SUBTRACT, Ld.getDeclaredIDs().get(0).clone())));
                Statement stmntRd = new ExpressionStatement(new AssignmentExpression(Rd.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, new BinaryExpression(Id.getDeclaredIDs().get(0).clone(), BinaryOperator.MODULUS, new NameID("size"))));
                badContainer.addStatement(stmntId);
                badContainer.addStatement(stmntRd);
                //b compute data partition bounds [CD_L_p,CD_U_p) owned by each p
        //if(p==0){CD_L_p=Ld;CD_U_p=CD_L_p+Id/size+Rd}else{CD_L_p=Ld+rank*Id/size+Rd);CD_U_p=CD_L_p+Ic/size}
        IfStatement cdpb = new IfStatement(new BinaryExpression(p.getDeclaredIDs().get(0).clone(), BinaryOperator.COMPARE_EQ, new IntegerLiteral(0)), new CompoundStatement());
        CompoundStatement cdpb_true_clause = new CompoundStatement();
        CompoundStatement cdpb_false_clause = new CompoundStatement();
        cdpb.setThenStatement(cdpb_true_clause);
        cdpb.setElseStatement(cdpb_false_clause);
        //for the true clause - master node
        Statement stmntCD_L_p_0 = new ExpressionStatement(new AssignmentExpression(CD_L_p.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, Ld.getDeclaredIDs().get(0).clone()));
            //Id_size = 'Id/size'
        Expression Id_size = new BinaryExpression(Id.getDeclaredIDs().get(0).clone(), BinaryOperator.DIVIDE, new NameID("size"));
        Statement stmntCD_U_p_0 = new ExpressionStatement(new AssignmentExpression(CD_U_p.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, new BinaryExpression(new BinaryExpression(CD_L_p.getDeclaredIDs().get(0).clone(),BinaryOperator.ADD,Id_size.clone()), BinaryOperator.ADD,Rd.getDeclaredIDs().get(0).clone() )));
        cdpb_true_clause.addStatement(stmntCD_L_p_0);
        cdpb_true_clause.addStatement(stmntCD_U_p_0);
        //for the false clause - other than master node ; note that p- indicates - rank - or index each processors
        //rank_Id_size = p*Id/size
        Expression rank_Id_size = new BinaryExpression(p.getDeclaredIDs().get(0).clone(), BinaryOperator.MULTIPLY, Id_size.clone());
        Statement stmntCD_L_p_1 = new ExpressionStatement(new AssignmentExpression(CD_L_p.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL,new BinaryExpression(new BinaryExpression(Ld.getDeclaredIDs().get(0).clone(),BinaryOperator.ADD,rank_Id_size.clone()), BinaryOperator.ADD,Rd.getDeclaredIDs().get(0).clone())));
        Statement stmntCD_U_p_1 = new ExpressionStatement(new AssignmentExpression(CD_U_p.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL,new BinaryExpression(CD_L_p.getDeclaredIDs().get(0).clone(), BinaryOperator.ADD, Id_size.clone())));
        cdpb_false_clause.addStatement(stmntCD_L_p_1);
        cdpb_false_clause.addStatement(stmntCD_U_p_1);
        badContainer.addStatement(cdpb);
        //c- get size of partition of CD owned by p
        //PSM = ..getPartitionSizeMultiplier
        Statement stmntPSM = new ExpressionStatement(new AssignmentExpression(PSM.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, broadcastADL1.getPartitionSizeMultiplier().clone()));
        //CD_S_p = (CD_U_p-CD_L_p)*PSM
        Statement stmntCD_S_p = new ExpressionStatement(new AssignmentExpression(CD_S_p.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, new BinaryExpression(new BinaryExpression(CD_U_p.getDeclaredIDs().get(0).clone(), BinaryOperator.SUBTRACT, CD_L_p.getDeclaredIDs().get(0).clone()), BinaryOperator.MULTIPLY, PSM.getDeclaredIDs().get(0).clone())));
        badContainer.addStatement(stmntPSM);
        badContainer.addStatement(stmntCD_S_p);
                
                //create MPI_Bcast(&AD[CD_L_p],CD_S_p , datatype, p, MPI_COMM_WORLD); statement
                FunctionCall MPI_Bcast = new FunctionCall(SymbolTools.getOrphanID("MPI_Bcast"));
                //add the arguments       
                 String argument1 = "&" + broadcastADL1.getVarExprWithParLpIndexReplaced(CD_L_p.getDeclaredIDs().get(0).clone()).clone();
                 MPI_Bcast.addArgument(new NameID(argument1));
                 MPI_Bcast.addArgument(CD_S_p.getDeclaredIDs().get(0).clone());
                //third argument 
                 String dataType = "MPI_" + broadcastADL1.getDataType().toString().toUpperCase();
                 MPI_Bcast.addArgument(new NameID(dataType));
                 MPI_Bcast.addArgument(p.getDeclaredIDs().get(0).clone());
                 MPI_Bcast.addArgument(new NameID("MPI_COMM_WORLD"));
                 Statement MPI_BcastStmnt = new ExpressionStatement(MPI_Bcast);
                 //add MPI fun to the container
                 badContainer.addStatement(MPI_BcastStmnt);
                 //add the container to the forloop body
                 broadcastLoopBody.addStatement(badContainer);
            }
            
            //add the forloop to the container
            container.addStatement(broadcastLoop);
        }//end of broadcastADL -comm code
        //step 4 comm code for broadcastADL
        if(!gatherADL.isEmpty()){
            //iterate throgh each gatherADL
            TranslationUnit tu = MPI_Utility.getTranslationUnit(program);
            for (MPI_AccessedData gatherADL1 : gatherADL){
                //declare rcvBuff for the gatherADL1 - global variable
                String varName = "rcvBuff"+gatherADL1.getVariableSymbol().getSymbolName();  
                VariableDeclaration rcvBuff = new VariableDeclaration(gatherADL1.getDataType(), new VariableDeclarator(new NameID(varName), gatherADL1.getVariableSymbol().getArraySpecifiers()));
                rcvBuff.annotate(new CommentAnnotation("rcv buffer for accessed data - "+gatherADL1.getVariableSymbol().getSymbolName()));
                //add the variables before the main procedue 
                tu.addDeclarationBefore(MPI_Utility.getMainProcedure(program),rcvBuff);
                
                //creat container - for the gatherADL1 - comms
                CompoundStatement gadContainer = new CompoundStatement();
                gadContainer.annotate(new CommentAnnotation("For accessed data - "+gatherADL1.getVariableSymbol().getSymbolName()+"- gathering"));
                //step 4.1 add statments - that compute comm variables
                //first find the AD partitioned dimension rangeAccess lb & ub
                Expression lb = gatherADL1.getAccessRange(gatherADL1.getDimOfIndexOftheParallelLoop()).getLB().clone();
                Expression ub = gatherADL1.getAccessRange(gatherADL1.getDimOfIndexOftheParallelLoop()).getUB().clone();
                Statement stmntLd = new ExpressionStatement(new AssignmentExpression(Ld.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, lb.clone()));
                Statement stmntUd = new ExpressionStatement(new AssignmentExpression(Ud.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, ub.clone()));
                gadContainer.addStatement(stmntLd);
                gadContainer.addStatement(stmntUd);
                //a- compute size of coomunication data and data partition remains
                //Id = Ud-Ld; Rd = Id%size;
                Statement stmntId = new ExpressionStatement(new AssignmentExpression(Id.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, new BinaryExpression(Ud.getDeclaredIDs().get(0).clone(), BinaryOperator.SUBTRACT, Ld.getDeclaredIDs().get(0).clone())));
                Statement stmntRd = new ExpressionStatement(new AssignmentExpression(Rd.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, new BinaryExpression(Id.getDeclaredIDs().get(0).clone(), BinaryOperator.MODULUS, new NameID("size"))));
                gadContainer.addStatement(stmntId);
                gadContainer.addStatement(stmntRd);
                //get partition size multiplier
                //PSM = ..getPartitionSizeMultiplier
                Statement stmntPSM = new ExpressionStatement(new AssignmentExpression(PSM.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, gatherADL1.getPartitionSizeMultiplier().clone()));
                gadContainer.addStatement(stmntPSM);
                //b compute data partition bounds [CD_L_p,CD_U_p) owned by each p
        //if(rank==0){CD_L_p=Ld;CD_U_p=CD_L_p+Id/size+Rd}else{CD_L_p=Ld+rank*Id/size+Rd);CD_U_p=CD_L_p+Ic/size}
        IfStatement cdpb = new IfStatement(new BinaryExpression(new NameID("rank"), BinaryOperator.COMPARE_EQ, new IntegerLiteral(0)), new CompoundStatement());
        CompoundStatement cdpb_true_clause = new CompoundStatement();
        CompoundStatement cdpb_false_clause = new CompoundStatement();
        cdpb.setThenStatement(cdpb_true_clause);
        cdpb.setElseStatement(cdpb_false_clause);
        //for the true clause - master node
        Statement stmntCD_L_p_0 = new ExpressionStatement(new AssignmentExpression(CD_L_p.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, Ld.getDeclaredIDs().get(0).clone()));
            //Id_size = 'Id/size'
        Expression Id_size = new BinaryExpression(Id.getDeclaredIDs().get(0).clone(), BinaryOperator.DIVIDE, new NameID("size"));
        Statement stmntCD_U_p_0 = new ExpressionStatement(new AssignmentExpression(CD_U_p.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, new BinaryExpression(new BinaryExpression(CD_L_p.getDeclaredIDs().get(0).clone(),BinaryOperator.ADD,Id_size.clone()), BinaryOperator.ADD,Rd.getDeclaredIDs().get(0).clone() )));
        cdpb_true_clause.addStatement(stmntCD_L_p_0);
        cdpb_true_clause.addStatement(stmntCD_U_p_0);
        //for the false clause - other than master node
        //rank_Id_size = rank*Id/size
        Expression rank_Id_size = new BinaryExpression(new NameID("rank"), BinaryOperator.MULTIPLY, Id_size.clone());
        Statement stmntCD_L_p_1 = new ExpressionStatement(new AssignmentExpression(CD_L_p.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL,new BinaryExpression(new BinaryExpression(Ld.getDeclaredIDs().get(0).clone(),BinaryOperator.ADD,rank_Id_size.clone()), BinaryOperator.ADD,Rd.getDeclaredIDs().get(0).clone())));
        Statement stmntCD_U_p_1 = new ExpressionStatement(new AssignmentExpression(CD_U_p.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL,new BinaryExpression(CD_L_p.getDeclaredIDs().get(0).clone(), BinaryOperator.ADD, Id_size.clone())));
        cdpb_false_clause.addStatement(stmntCD_L_p_1);
        cdpb_false_clause.addStatement(stmntCD_U_p_1);
        gadContainer.addStatement(cdpb);
                
                //if Rd is !=0, copy gatherADL1 to the buffer, in the master processor
                //prepare if(rank==0 && Rd!=0){} -- partRmDataCopy
                Expression partRmDataCopy_condition = new BinaryExpression(new BinaryExpression(new NameID("rank"), BinaryOperator.COMPARE_EQ, new IntegerLiteral(0)),BinaryOperator.LOGICAL_AND,new BinaryExpression(Rd.getDeclaredIDs().get(0).clone(), BinaryOperator.COMPARE_NE, new IntegerLiteral(0)));
                CompoundStatement partRmDataCopy_body = new CompoundStatement();
                IfStatement partRmDataCopy_if = new IfStatement(partRmDataCopy_condition, partRmDataCopy_body);
                      //update CD_L_p - parti starting point  
                    //prepare  'CD_L_p=CD_L_p+Rd);' statement - because - master node doesnt include partRmData in the gathering comm
                //Expression Rd_PSM = new BinaryExpression(Rd.getDeclaredIDs().get(0).clone(), BinaryOperator.MULTIPLY,PSM.getDeclaredIDs().get(0).clone());
                Statement stmntCD_L_p_update = new ExpressionStatement(new AssignmentExpression(CD_L_p.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, new BinaryExpression(CD_L_p.getDeclaredIDs().get(0).clone(), BinaryOperator.ADD,Rd.getDeclaredIDs().get(0).clone())));
                partRmDataCopy_body.addStatement(stmntCD_L_p_update);
		//prepare memcpy function satement-to copy the accesseddata content to the buffer
		//memcpy(rcvBuffAD,AD,sizeof(AD));
                FunctionCall memcpy_fc = new FunctionCall(new NameID("memcpy"));
                memcpy_fc.addArgument(rcvBuff.getDeclaredIDs().get(0).clone());
                memcpy_fc.addArgument(new NameID(gatherADL1.getVariableSymbol().getSymbolName()));
                //prepare the third argument - sizeof(AD)
                FunctionCall sizeof_fc = new FunctionCall(new NameID("sizeof"));
                sizeof_fc.addArgument(new NameID(gatherADL1.getVariableSymbol().getSymbolName()));
                memcpy_fc.addArgument(sizeof_fc);
                Statement memcpy_stmnt = new ExpressionStatement(memcpy_fc);
                partRmDataCopy_body.addStatement(memcpy_stmnt);
                //add partRmDataCopy_if to the gadcontainer
                gadContainer.addStatement(partRmDataCopy_if);
                //c- get size of partition of CD owned by p
         //CD_S_p = (CD_U_p-CD_L_p)*PSM
        Statement stmntCD_S_p = new ExpressionStatement(new AssignmentExpression(CD_S_p.getDeclaredIDs().get(0).clone(), AssignmentOperator.NORMAL, new BinaryExpression(new BinaryExpression(CD_U_p.getDeclaredIDs().get(0).clone(), BinaryOperator.SUBTRACT, CD_L_p.getDeclaredIDs().get(0).clone()), BinaryOperator.MULTIPLY, PSM.getDeclaredIDs().get(0).clone())));
        gadContainer.addStatement(stmntCD_S_p);
                //prepare the MPI gather function
                //gather to the buffer
		//MPI_Gather(&AD[CD_L_p], CD_S_p, Datatype,&rcvBuffAD[CD_L_p], CD_S_p,Datatype, 0,MPI_COMM_WORLD);
                FunctionCall MPI_Gather = new FunctionCall(SymbolTools.getOrphanID("MPI_Gather"));
                //add the arguments       
                 String argument1 = "&" + gatherADL1.getVarExprWithParLpIndexReplaced(CD_L_p.getDeclaredIDs().get(0).clone()).clone();
                 MPI_Gather.addArgument(new NameID(argument1));
                 MPI_Gather.addArgument(CD_S_p.getDeclaredIDs().get(0).clone());
                //third argument 
                 String dataType = "MPI_" + gatherADL1.getDataType().toString().toUpperCase();
                 MPI_Gather.addArgument(new NameID(dataType));
                 //prepare the 4th argument - from the argument1
                 //simply add 'rcvBuff' before the AD
                 String argument4 = "&" +"rcvBuff"+ gatherADL1.getVarExprWithParLpIndexReplaced(CD_L_p.getDeclaredIDs().get(0).clone()).clone();
                 MPI_Gather.addArgument(new NameID(argument4));
                 MPI_Gather.addArgument(CD_S_p.getDeclaredIDs().get(0).clone());
                 //6th argument 
                 MPI_Gather.addArgument(new NameID(dataType));
                 //7th argu - destin
                 MPI_Gather.addArgument(new IntegerLiteral(0));
                 MPI_Gather.addArgument(new NameID("MPI_COMM_WORLD"));
                 Statement MPI_GatherStmnt = new ExpressionStatement(MPI_Gather);
                 //add MPI_GatherStmnt to the gadcontainer
                gadContainer.addStatement(MPI_GatherStmnt);
                
                //copy from the buffer to the AD - on the master node - buff2AD
                BinaryExpression buff2AD_condi = new BinaryExpression(new NameID("rank"), BinaryOperator.COMPARE_EQ, new IntegerLiteral(0));
                IfStatement buff2AD_if = new IfStatement(buff2AD_condi, new CompoundStatement());
                //prepare memcpy function satement-to copy the buffer data to accesseddata variable
		//memcpy(rcvBuffAD,AD,sizeof(AD));
                FunctionCall memcpy_fc_2 = new FunctionCall(new NameID("memcpy"));
                memcpy_fc_2.addArgument(new NameID(gatherADL1.getVariableSymbol().getSymbolName()));
                memcpy_fc_2.addArgument(rcvBuff.getDeclaredIDs().get(0).clone());
                //prepare the third argument - sizeof(rcvBuffAD)
                FunctionCall sizeof_fc_2 = new FunctionCall(new NameID("sizeof"));
                sizeof_fc_2.addArgument(rcvBuff.getDeclaredIDs().get(0).clone());
                memcpy_fc_2.addArgument(sizeof_fc_2);
                Statement memcpy_stmnt_2 = new ExpressionStatement(memcpy_fc_2);
                //add memcpy_stmnt_2 to the buff2AD_if body
                buff2AD_if.setThenStatement(memcpy_stmnt_2);
                //add buff2AD_if to the gadcontainer
                gadContainer.addStatement(buff2AD_if);
                
                //add the gadContainer to the main container
                container.addStatement(gadContainer);
            }
        }//end of gatherADL -comm code
       
        /////add time instrumenting at the end of comm////
        if(timeReport){
            //prepare mpipar_comm_timer_stop(); statement
            FunctionCall stop_timer=  new FunctionCall(SymbolTools.getOrphanID("mpipar_comm_timer_stop"));
            Statement stop_timer__stmnt = new ExpressionStatement(stop_timer);
            container.addStatement(stop_timer__stmnt);
        }
       
        
        //////////////////////
        //finally add the comm ifcontainer after the loop
        CompoundStatement parent = (CompoundStatement) loop.getParent();
        parent.addStatementAfter(loop, ifContainer);
    }
    
    /**
     * Insert time instrumenting functions
     * @param program 
     */
    public void insertTimeInstrumentation(Program program){
        //Find #pragma scop and #pragma endscop annotations --to spot the kerenl - in the program main function
        //iterte through each line in the main proc body and get 
        Statement kernelBegin=null,kernelEnd=null;
        FlatIterator flt_itr_annot = new FlatIterator(MPI_Utility.getMainProcedure(program).getBody());
        
        if(flt_itr_annot!=null){
            while (flt_itr_annot.hasNext()) {
                Object next = flt_itr_annot.next();
                //System.out.println("flat iterator object--"+next);
                if(next instanceof Statement){
                    Statement annot = (Statement) next;
                    if(annot.toString().contains("#pragma scop"))
                        kernelBegin = annot;
                     if(annot.toString().contains("#pragma endscop"))
                        kernelEnd = annot;
                   
                }
            }
        }
        
        //after spoting the kernel begin and end point --insert timing functions
        if(kernelBegin!=null && kernelEnd!=null){
            //prepare polybench_timer_start(); function call
            FunctionCall start_timer=  new FunctionCall(SymbolTools.getOrphanID("polybench_timer_start"));
            Statement time_start_stmnt = new ExpressionStatement(start_timer);
            //insert before kernel begins
            MPI_Utility.getMainProcedure(program).getBody().addStatementBefore(kernelBegin, time_start_stmnt);
            //prepare  polybench_timer_stop(); function call  
            FunctionCall stop_timer=  new FunctionCall(SymbolTools.getOrphanID("polybench_timer_stop"));
            Statement time_stop_stmnt = new ExpressionStatement(stop_timer);
            //insert before kernel end statement
            MPI_Utility.getMainProcedure(program).getBody().addStatementBefore(kernelEnd, time_stop_stmnt);
            //prepare if(rank==0) mpipar_comm_overhead_print(); statement
            //prepare rank==0 expression
            Expression condition = new BinaryExpression(new NameID("rank"), BinaryOperator.COMPARE_EQ, new IntegerLiteral(0));
            IfStatement if_rank_0 = new IfStatement(condition, new CompoundStatement());
            //prepare  mpipar_comm_overhead_print(); function call
            FunctionCall print_comm_overhead=  new FunctionCall(SymbolTools.getOrphanID("mpipar_comm_overhead_print"));
            Statement print_comm_overhead_stmnt = new ExpressionStatement(print_comm_overhead);
            if_rank_0.setThenStatement(print_comm_overhead_stmnt);
            //insert after kernel end statement
            MPI_Utility.getMainProcedure(program).getBody().addStatementAfter(kernelEnd, if_rank_0);
        }
        
       
    }

}
