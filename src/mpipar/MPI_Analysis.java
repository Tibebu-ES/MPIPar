/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mpipar;

import cetus.analysis.AnalysisPass;
import cetus.analysis.LoopTools;
import cetus.analysis.RangeAnalysis;
import cetus.analysis.RangeDomain;
import cetus.hir.Annotation;
import cetus.hir.AssignmentExpression;
import cetus.hir.AssignmentOperator;
import cetus.hir.BinaryExpression;
import cetus.hir.BinaryOperator;
import cetus.hir.CetusAnnotation;
import cetus.hir.CommentAnnotation;
import cetus.hir.DepthFirstIterator;
import cetus.hir.Expression;
import cetus.hir.ExpressionStatement;
import cetus.hir.ForLoop;
import cetus.hir.PragmaAnnotation;
import cetus.hir.Program;
import cetus.hir.Specifier;
import cetus.hir.Statement;
import cetus.hir.Symbol;
import cetus.hir.SymbolTools;
import cetus.hir.Traversable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

/**
 * find mpi parallel loops and shared variables info with in these loops attach
 * mpi parallel and mpi sharedvariables (list<SharedVariables >) annotations to
 * found parallel loops also add mpi reduction (map<String op, Set<Expression>>)
 * annotation to the loops whic are used by the MPI_LoopTransform pass
 *
 * @author tibebu
 */
public class MPI_Analysis extends AnalysisPass {

    public MPI_Analysis(Program program) {
        super(program);
    }

    @Override
    public String getPassName() {
        return "MPI Analysis";
    }

    @Override
    public void start() {
        getMPIParallelLoops(program);
    }

    /**
     * look in main procedure parse internal annotations and find cetus parallel
     * loops & MPI normalized loops find shared variables with in this parallel
     * loop add MPI parallel annotation to the parallel loops add MPI
     * accessed_data(list of variables) annotation to the parallel loop
     *
     * @param program
     */
    public void getMPIParallelLoops(Program program) {
        int numberOfParallelLoopsFound = 0;
        int numberOfOuterParallelLoopsFound = 0;
        //for tracking the last outer paralle loop - this loop doesnt need second comm code -hence #MPI nsc will be annoted to it
        Stack<ForLoop> parallelLoopsStack = new Stack<>();
        DepthFirstIterator bfs_itr_lp = new DepthFirstIterator(MPI_Utility.getMainProcedure(program));
  //      DepthFirstIterator bfs_itr_lp = null;
   //      if(MPI_Utility.getKernelProcedure(program)!=null)
    //     bfs_itr_lp = new DepthFirstIterator(MPI_Utility.getKernelProcedure(program));
        if (bfs_itr_lp != null) {
            while (bfs_itr_lp.hasNext()) {
                //get the 
                Object o = bfs_itr_lp.next();
                if (o instanceof ForLoop) {
                    ForLoop loop = (ForLoop) o;
                    //find cetus parallel & mpi normalized loops
                    if (loop.containsAnnotation(CetusAnnotation.class, "parallel") && loop.containsAnnotation(MPI_Annotation.class, MPI_Annotation.KEY_NORMALIZED)) {
          //           if (loop.containsAnnotation(CetusAnnotation.class, "parallel")) {                
                        //get the loop depth - if outer or inner -- by parsing the loop name and
                        //counting # in the name- if 
                        //Annotation loopName = lo
                        //get the loop depth
                        String _lname = "";
                        int loopDepth = 0;
                        if (loop.containsAnnotation(PragmaAnnotation.class, "name")) {
                            PragmaAnnotation loopann = loop.getAnnotation(PragmaAnnotation.class, "name");
                            String lname = loopann.get("name");
                            _lname = lname;
                            String[] tokens = lname.split("#");
                            int numOfHashTag = tokens.length - 1;
                            loopDepth = numOfHashTag - 1;

                        }

                        numberOfParallelLoopsFound++;
                        if (loopDepth == 0) {
                            numberOfOuterParallelLoopsFound++;
                        }
                        //attache MPI parallel annotation to the target loop
                        MPI_Annotation mpiAnn = new MPI_Annotation(MPI_Annotation.KEY_PARALLEL, loopDepth);
                        loop.annotate(mpiAnn);

                        //put the reduction variables info - create MPI reduction(MAP<String operator, Set<Expression>>)
                        if (loop.containsAnnotation(CetusAnnotation.class, "reduction")) {
                            CetusAnnotation cetusRedAnno = loop.getAnnotation(CetusAnnotation.class, "reduction");
                            loop.annotate(new MPI_Annotation(MPI_Annotation.KEY_REDUCTION, cetusRedAnno.get("reduction")));
                        }

                        //compute accessed data
                       
                        List<MPI_AccessedData> accessedVariables = findAccessedDataInAparallelLoop(loop);

                        //create MPI accessedDataLst() annotation
                        if (!accessedVariables.isEmpty()) {
                            MPI_Annotation mpiSvAnno = new MPI_Annotation(MPI_Annotation.KEY_ACCESSED_DATA, accessedVariables);
                            loop.annotate(mpiSvAnno);
                        }
                        
                        //push the loop to the stack - if it is outer loop
                        if(loopDepth==0)
                            parallelLoopsStack.push(loop);
                        
                    }

                }
            }
        }
        //get  array of outerparallel loops - in lexical order , i.e first element comes first lexically
         ForLoop outerParallelLoops[] = new ForLoop[numberOfOuterParallelLoopsFound];
         for(int i=numberOfOuterParallelLoopsFound-1;i>=0;i--){
             outerParallelLoops[i]=parallelLoopsStack.pop();
         }
        //determine type of comm required by each accessed data in each parallel loop
         setTypeOfCommunication(outerParallelLoops);
         //add analysis summary- comment at the end of the source code
        //get the benchmark name - from the kernel proc name which is 'kernel_xxx' where xxx=benchname
        String benchName = "";
        if (MPI_Utility.getKernelProcedure(program) != null) {
            String kernelProcName = MPI_Utility.getKernelProcedure(program).getName().toString();
            String tokens[] = kernelProcName.split("_");
            if (tokens.length > 1) {
                benchName = tokens[1];
            }
        }

        Annotation summary = new CommentAnnotation("=>========================MPI Analysis===================\n"
                + "=> Benchmark =" + benchName + "\n"
                + "=> TOTAL NO. OF PARALLEL LOOPS FOUND = " + numberOfParallelLoopsFound + "\n"
                + "=> INNER PARALLEL LOOPS FOUND =" + (numberOfParallelLoopsFound - numberOfOuterParallelLoopsFound) + "\n"
                + "=> OUTER PARALLEL LOOPS FOUND = " + numberOfOuterParallelLoopsFound + "\n"
                + "=>==========================================================");
        MPI_Utility.getMainProcedure(program).annotateAfter(summary);

    }

    /**
     * find list of accessed data with in the given parallel floop given the
     * forloop iterate each child of the floop body, check if itis Assignment
     * expression- get the leftSide variable chck if it array
     *
     * floop index can be anywhere - but locate where it is get the variable
     * datataype- create new sharedVariable object whether the variable is
     * reduction variable or not make sure no duplicate
     *
     * @param floop
     *
     * @return
     */
    public List<MPI_AccessedData> findAccessedDataInAparallelLoop(ForLoop floop) {
        List<MPI_AccessedData> accessedDataLst = new ArrayList<>();
        Map<Statement, RangeDomain> range_dom = RangeAnalysis.getRanges(floop);
        Statement loopBody = floop.getBody();
        Symbol index = LoopTools.getLoopIndexSymbol(floop);
        DepthFirstIterator childrens = new DepthFirstIterator(loopBody);
        while (childrens.hasNext()) {
            Traversable next = childrens.next();
            //check it is assignment Statment
            if (next instanceof ExpressionStatement) {
                ExpressionStatement exprSta = (ExpressionStatement) next;
                if (exprSta.getExpression() instanceof AssignmentExpression) {
                    AssignmentExpression assExpr = (AssignmentExpression) exprSta.getExpression();
                    MPI_AccessedData accssedData = new MPI_AccessedData(SymbolTools.getSymbolOf(assExpr.getLHS()));
                     //set the parallelil parent floop index
                    accssedData.setIndexOftheParallelLoop(index);
                    //check if it is reduction variable or not -- by parsing the cetus annotaion
                    if (floop.containsAnnotation(CetusAnnotation.class, "reduction")) {
                        System.out.println("reduction found------------------------------>");
                        CetusAnnotation cetusReduAnnt = floop.getAnnotation(CetusAnnotation.class, "reduction");
                        Map<String, Set<Expression>> redMap = cetusReduAnnt.get("reduction");
                        //System.out.println("reduction found"+redMap.toString());
                        for (String op : redMap.keySet()) {
                            Set<Expression> redExpres = redMap.get(op);
                            if (redExpres.contains(assExpr.getLHS())) {
                                accssedData.setReductionOperator(BinaryOperator.fromString(op));
                            }
                        }
                    }//checking for reduction variable - setting  finished 
                    
                    //check if accssed data is scalar or array
                    if (SymbolTools.isScalar(accssedData.getVariableSymbol())) {
                        //if it is scalar  -- fill infos related to scalar 
                        accssedData.setIsArray(false);
                        //no need to set the var expression & range domain                    
                        //check if it is reduction variable -if not then no need to share this among process
                        //added to solve - the problem that occur in parallel floop that has inner floop
                        //in this case the pass  consider inner loops index as a shared variable - which shouldnt be
                        //hence i add this restriction that- if it is not reduction --then it is not shared variable
                        if (accssedData.isIsReductionVariable()) {
                            accessedDataLst.add(accssedData);
                        }
                    } else if (SymbolTools.isArray(accssedData.getVariableSymbol())) {
                        //if the variable is array
                        accssedData.setIsArray(true);
                        accssedData.setVarExpression(assExpr.getLHS().clone());
                        //set the rangeDomain - holds the range for all symbols in the varExpression
                        accssedData.setRangeDomain(range_dom.get(assExpr.getLHS().getStatement()));
                        //check if accssed data is already in the list - if not add it to the list
                        if (!accessedDataLst.contains(accssedData)) 
                            accessedDataLst.add(accssedData);
                    } //if  array - finished
                }//end of if Ass expre
            }//end of if expre statment
        }

        return accessedDataLst;
    }
    
    /**
     * determine type of comm needed by each accessed data in each outer parallel loop
     * gathering or broadcasting 
     * @param ploops 
     */
    public void setTypeOfCommunication(ForLoop oploops[]){
        //iterate throgh each loop
        for(int i=0;i<oploops.length;i++){
            //System.out.println("ACCESSED SYMBOLS-"+SymbolTools.getAccessedSymbols(oploops[i]));
            //get list of accessed data
            if(oploops[i].containsAnnotation(MPI_Annotation.class, MPI_Annotation.KEY_ACCESSED_DATA)){
                MPI_Annotation mpiSVannot = oploops[i].getAnnotation(MPI_Annotation.class, MPI_Annotation.KEY_ACCESSED_DATA);
                List<MPI_AccessedData> accessedVariables = mpiSVannot.get(MPI_Annotation.KEY_ACCESSED_DATA);
                //iterate through each accesseddata
                for (MPI_AccessedData accessedVariable : accessedVariables) {
                    //checks if accessedVariable is found in one of the prceding parallel loop
                    boolean notAccessedAnyWhere=true;
                    for(int j=i+1;j<oploops.length;j++){
                        //get accessed symbols
                        Set<Symbol> accessedSymbols =SymbolTools.getAccessedSymbols(oploops[j]);
                        
                        if(accessedSymbols.contains(accessedVariable.getVariableSymbol())){
                            notAccessedAnyWhere=false;
                            break;
                        }
                    }
                    //set the commTyep - 
                    if(notAccessedAnyWhere)
                        accessedVariable.setCommType(MPI_AccessedData.COMMUNICATION_TYPE_GATHERING);
                    else
                        accessedVariable.setCommType(MPI_AccessedData.COMMUNICATION_TYPE_BROADCASTING);
                    
                    //for debugging -
                   // System.out.println("accessedData-->"+accessedVariable);
                }

            }
            
        }
    }

}
