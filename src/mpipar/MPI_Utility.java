/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mpipar;

import cetus.analysis.DDGraph;
import cetus.analysis.LoopTools;
import cetus.analysis.RangeAnalysis;
import cetus.analysis.RangeDomain;
import cetus.hir.Annotatable;
import cetus.hir.Annotation;
import cetus.hir.AnnotationDeclaration;
import cetus.hir.AssignmentExpression;
import cetus.hir.AssignmentOperator;
import cetus.hir.BinaryExpression;
import cetus.hir.BinaryOperator;
import cetus.hir.BreadthFirstIterator;
import cetus.hir.CetusAnnotation;
import cetus.hir.CommentAnnotation;
import cetus.hir.CompoundStatement;
import cetus.hir.Declaration;
import cetus.hir.Declarator;
import cetus.hir.DepthFirstIterator;
import cetus.hir.Expression;
import cetus.hir.ExpressionStatement;
import cetus.hir.ForLoop;
import cetus.hir.FunctionCall;
import cetus.hir.IDExpression;
import cetus.hir.Identifier;
import cetus.hir.IntegerLiteral;
import cetus.hir.Loop;
import cetus.hir.NameID;
import cetus.hir.Procedure;
import cetus.hir.ProcedureDeclarator;
import cetus.hir.Program;
import cetus.hir.Specifier;
import cetus.hir.Statement;
import cetus.hir.StatementExpression;
import cetus.hir.Symbol;
import cetus.hir.SymbolTable;
import cetus.hir.SymbolTools;
import cetus.hir.Symbolic;
import cetus.hir.TranslationUnit;
import cetus.hir.Traversable;
import cetus.hir.UnaryExpression;
import cetus.hir.UnaryOperator;
import cetus.hir.VariableDeclaration;
import cetus.hir.VariableDeclarator;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 *
 * @author tibebu
 */
public class MPI_Utility {

    public static void InstrumentLoops(Program p) {
        /*Create cetus_tic and cetus_toc identifiers*/
       
        Declarator d = new VariableDeclarator(new NameID("cetus_tic"));
        Identifier cetus_tic = SymbolTools.getOrphanID("cetus_tic");
        Identifier cetus_toc = SymbolTools.getOrphanID("cetus_toc");
        /*for each proc in p*/
        DepthFirstIterator dfs_itr = new DepthFirstIterator(p);
        while (dfs_itr.hasNext()) {
            Object op = dfs_itr.next();
            if (op instanceof Procedure) {
                Procedure proc = (Procedure) op;

                List<Statement> loops = new ArrayList<Statement>();
                /*for each loop in proc find loops*/
                DepthFirstIterator dfs_itr_lp = new DepthFirstIterator(proc);
                while (dfs_itr_lp.hasNext()) {
                    Object o = dfs_itr_lp.next();
                    if (o instanceof ForLoop) {
                        loops.add((Statement) o);
                    }
                } //Finding loop in proc finishes
                String name = proc.getSymbolName();
                int i = 0;
                for (Statement loop : loops) {
                    /*create annotation*/
                    Annotation loopName = new CommentAnnotation(name + "#" + i);
                    /* creates a timing call function*/
                    Statement ticstat = CreateTimingCall(cetus_tic, i);
                    Statement tocstat = CreateTimingCall(cetus_toc, i++);
                    AddInstrumentation(loop, loopName, ticstat, tocstat);
                }
            }
        }
    }

    /*Create timing funcation call 'funcName(num);'*/
    public static Statement CreateTimingCall(Identifier fcnName, int num) {
        FunctionCall fcnCall = new FunctionCall(fcnName.clone());
        fcnCall.addArgument(new IntegerLiteral(num));
        Statement fcnStat = new ExpressionStatement(fcnCall);
        return fcnStat;
    }
    /*Add annotation and tic statment before the loop, add toc statment after the loop*/

    public static void AddInstrumentation(Statement loop, Annotation loopname, Statement tic, Statement toc) {
        /*locate parent compound statment*/
        CompoundStatement parent = (CompoundStatement) loop.getParent();
        loop.annotate(loopname);
        parent.addStatementBefore(loop, tic);
        parent.addStatementAfter(loop, toc);
    }

    
   
    /**
     *  Inserts a variable at the given scope - t
     * variableSpecifier = Specifier.INT, 
     * @param t
     * @param variableName
     * @param variableSpecifier
     * @return return the declared variable
     */
    public static VariableDeclaration DeclareVariable(Traversable t, String variableName, Specifier variableSpecifier) {

        
        //find scope
        while (!(t instanceof SymbolTable)) {
            t = t.getParent();
        }
        SymbolTable symtab = (SymbolTable) t;
        //check the name validation
        NameID nameId = new NameID(variableName);
        int i = 0;
        while (symtab.findSymbol(nameId) != null) {
            nameId = new NameID(variableName + (i++));
        }
        //create declaration
        Declarator d = new VariableDeclarator(nameId);
        VariableDeclaration decl = new VariableDeclaration(variableSpecifier, d);
        //insert declaration at the end
        
        symtab.addDeclaration(decl);
        //return the declaration
        return decl;
       // System.out.println("Done!!!!!!!!!!!!!!!!!!!!!");
    }

    /**
     * return main procedure 
     * @param p
     * @return 
     */
    public static Procedure getMainProcedure(Program p) {
        Procedure main = null;
        //find the main procedure
        
        DepthFirstIterator dfs_itr = new DepthFirstIterator(p);
        while (dfs_itr.hasNext()) {
            Object op = dfs_itr.next();
            if (op instanceof Procedure) {
                Procedure proc = (Procedure) op;
                if (proc.getSymbolName().equalsIgnoreCase("main")) {
                    main = proc;
                }

            }
        }

        return main;
    }
    
    /**
     * find and return procedure starts with kernel
     * if not found return null
     * @param p
     * @return 
     */
    public static Procedure getKernelProcedure(Program p){
        Procedure kernel = null;
        BreadthFirstIterator dfs_itr = new BreadthFirstIterator(p);
        while (dfs_itr.hasNext()) {
            Object op = dfs_itr.next();
            if (op instanceof Procedure) {
                Procedure proc = (Procedure) op;
                if (proc.getSymbolName().startsWith("kernel")) {
                    kernel = proc;
                }

            }
        }
        return kernel;
    }
    
    /**
     * Get the program trnslationUnit
     * @param p
     * @return 
     */
    public static TranslationUnit getTranslationUnit(Program p){
        //get the program trnslation unit
        TranslationUnit tu = null;
        DepthFirstIterator dfs_itr = new DepthFirstIterator(p);
        while(dfs_itr.hasNext()){
            Object o = dfs_itr.next();
            if(o instanceof TranslationUnit){
                tu = (TranslationUnit)o;
                break;
            }
        }
        return tu;
    }

    /**
     * search expression a with in t and replace them with b
     * @param a - expression to be replaced - old expression
     * @param b - new expression
     * @param t 
     */
    public static void FindAndReplaceExpressions(Expression a,Expression b, Traversable t){
        List oldExprs = new ArrayList();
        DepthFirstIterator dpfit = new DepthFirstIterator(t);
        //search
        while (dpfit.hasNext()) {
            Object child = dpfit.next();
            if(child.equals(a))
                oldExprs.add(child);
        }
        //replace
        for (Object oldExpr : oldExprs) {
            Expression newExp = b.clone();
            ((Expression)oldExpr).swapWith(newExp);
        }
    }
    
    /* using dependence info for parallelization*/
    /**
     * check if a given loop can be parallelized
     * check if scalar dependence present using- privatization and reduction variable information
     * check if array carried dependence presents
     * @param pdg
     * @param loop
     * @return 
     */
    public static boolean IsLoopParallel(DDGraph pdg,Loop loop) {
        boolean isParallel = false;
        
        //check eligibility for dependence testing
        if (LoopTools.checkDataDependenceEligibility(loop) == true) {
            /*check if scalar depence present using privatization and reduction variable information*/
            if (LoopTools.scalarDependencePossible(loop) == true) {
                isParallel = false;
            } //check if array carried loop dependence exit
            else if (pdg.checkLoopCarriedDependence(loop) == true) {
                isParallel = false;
            } //no dependence
            else {
                isParallel = true;
            }
        }
        return isParallel;
    }

    /**
     * iterate through all the forloops in the program - try to normaliz it - 
     * annotate if it is normalized or not - above each for loop
     * @param p 
     */
    public static void DetermineNormalizedLoops(Program p){
        DepthFirstIterator it = new DepthFirstIterator(p);
        while (it.hasNext()) {
            Object next = it.next();
            if(next instanceof ForLoop){
                ForLoop nextFl = (ForLoop) next;
                boolean isNormalized = NormalizeLoop(nextFl);
                nextFl.annotate(new CommentAnnotation("Is Normalized ="+isNormalized));
            }
            
        }
    }
    /**
     * 
     * determine parallel loops and add annotation[whether loop is parallel or not] on each loops
     * @param program
     */
    public static void DetermineParallelLoops(Program program) {
           DDGraph dpg = program.getDDGraph();
           if(dpg==null){
               program.createNewDDGraph();
               dpg=program.getDDGraph();
           }
           //print ddgraph
         //  System.out.println("DDGraph->"+dpg.toString());
            /*for each loop in proc find loops*/
           DepthFirstIterator dfs_itr_lp = new DepthFirstIterator(program);
           while (dfs_itr_lp.hasNext()) {
               Object o = dfs_itr_lp.next();
               if (o instanceof ForLoop) {
                   Loop loop = (Loop)o;
                   //check if it is parallel
                   boolean isParallel = IsLoopParallel(dpg,loop);
                  // System.out.println("Is Parallel ="+isParallel);
                   //add annotation
                   Annotation info = new CommentAnnotation("Is Parallel ="+isParallel);
                   ((Statement)loop).annotate(info);
                   
                         
               }
           }   
       
        
    }
    
    
    /**
     * Iterate through all the forloops with in the main procedure - check if it can be parallelized 
     * return list of parallel forloops
     * also annotate each forloop -
     * @param program
     * @return 
     */
    public static List<Loop> GetParallelLoops(Program program){
        List<Loop> parallelLoops = new ArrayList<>();
        DDGraph dpg = program.getDDGraph();
        
           if(dpg==null){
               System.out.println("dpg was null");
               program.createNewDDGraph();
               dpg=program.getDDGraph();
           }
           
            /*for each loop in proc find loops*/
           DepthFirstIterator bfs_itr_lp = new DepthFirstIterator(getMainProcedure(program));
           while (bfs_itr_lp.hasNext()) {
               //get the 
               Object o = bfs_itr_lp.next();
               if (o instanceof ForLoop) {
                   Loop loop = (Loop)o;                  
                   //check if it is parallel
                   boolean isParallel = IsLoopParallel(dpg,loop);
                   if(isParallel){ 
                            parallelLoops.add(loop);
                   }
                   ((ForLoop)loop).annotate(new CommentAnnotation("Is Parallel"+isParallel));
               }
           }   

        
        return parallelLoops;
    }
    
    //NORMALIZE LOOP
    
    /**
     * try to normalize given loop and if succeeded return true if not false
     * normalizing the initial - result is ignored since it doesnt cause problem to the loop partitioning
     * @param loop
     * @return 
     */
    public static boolean NormalizeLoop(ForLoop loop){
         
        boolean status = true;
        //normalize the loop condition - initial and step
        boolean statusLC = normalizeLoopCondition(loop);
        boolean statusLI = normalizeLoopInitial(loop); // the result of this normalization is not important
        boolean statusLS = normalizeLoopStep(loop);
        status = statusLC && statusLS;
        //simplify 
        //Symbolic.simplifyWithin(loop);
        
        return status;
    }
    
    /**
     * first check if it is outerbounded -if not return false
     * if it is outerbounded - replace the outerbound expression with the outerbound
     * try to normalize given loop and if succeeded return true if not false
     * currently works on < & <= operators 
     * @param loop
     * @return 
     */
    public static boolean normalizeLoopCondition(ForLoop loop){
        
        boolean result = true;
        
        ForLoop normalizedLoop = loop;
        //1. get the loop condition expression
        Expression condition = normalizedLoop.getCondition();
        //2. check if the condition is binary expression and outerbounded- ifnot set the result false 
        if(LoopTools.isUpperBoundConstant(loop) && condition instanceof BinaryExpression){
            BinaryExpression bexpr = (BinaryExpression)condition;
            Expression ub = LoopTools.getUpperBoundExpression(loop).clone();
            //3 get the operator
            BinaryOperator op =bexpr.getOperator();
            if(op == BinaryOperator.COMPARE_LT || op == BinaryOperator.COMPARE_LE){
                //4 if the operator is < or <=  prepare the new loopconditions as loopIndex < (ub+1)
                //and set the loop condition 
                Expression newLpCondition = new BinaryExpression(LoopTools.getIndexVariable(loop).clone(),BinaryOperator.COMPARE_LT, Symbolic.simplify(new BinaryExpression(ub, BinaryOperator.ADD, new IntegerLiteral(1))));
                loop.setCondition(newLpCondition);
                result = true;
            }else{
                //if the operator is other than the above return false - which means cant be normalized or not normalized
                result = false;
            }
            
        }else{
            result = false;
        }
        return result;
    }
    
    /**
     * try to change the loop index inreament or step to 1
     * check if the loop increament is one - if so return true
     * if it is other than one try to normalize it to 1 - if succeede return true
     * if not succeded return false
     * @param loop
     * @return 
     */
    public static boolean normalizeLoopStep(ForLoop loop){
        boolean result = false;
        //get the loop index increament - step value
        int stepValue =getLoopStepValue(loop);
        if(stepValue == 1){
            //no need to structure the loop just - set the result trrue
            result = true;
        }else if(stepValue >1){
            //restrcture the loop
                //get the old upper bound
            Expression upbOld = Symbolic.simplify(FindLoopUpperBound(loop));
            //cretae the new upper bound as upbOld/stepValue
            //get the newupb integer value
            if(upbOld instanceof IntegerLiteral){
                int upbOldValue = (int)((IntegerLiteral)upbOld).getValue();
                int upbNewValue = upbOldValue/stepValue;
                //round to the next integer --if there is remainder
                if((upbOldValue%stepValue)>0)
                    upbNewValue+=1;
                
                //Expression upbNew = new BinaryExpression(upbOld.clone(), BinaryOperator.DIVIDE,new IntegerLiteral(stepValue));
                  Expression upbNew = new IntegerLiteral(upbNewValue);
                //simplify expression
                upbNew = Symbolic.simplify(upbNew);
                //set the new Upperbound
                setLoopUpperBound(loop, upbNew.clone());
                //prepare the stepExpres index++
                UnaryExpression stepEx = new UnaryExpression(UnaryOperator.POST_INCREMENT, FindLoopIndex(loop).clone());
                //set tthe step expression
                loop.setStep(stepEx.clone());
                
                //replace index exp i inside the loop by i*stepValue  - to reverse the change made above - 
                
                Expression oldExp = FindLoopIndex(loop);
                Expression newExp = new BinaryExpression(oldExp.clone(), BinaryOperator.MULTIPLY, new IntegerLiteral(stepValue));
                FindAndReplaceExpressions(oldExp, newExp, loop.getBody());
                /*
                
                //TODO update the oldindex value by adding the flg expression at the begining of the loop 
                //index = index*stepValue;
                Expression index = FindLoopIndex(loop);
                Expression updateIndexExp = new AssignmentExpression(index.clone(), AssignmentOperator.NORMAL, new BinaryExpression(index.clone(), BinaryOperator.MULTIPLY, new IntegerLiteral(stepValue)));
                Statement  updateIndexExpStat = new ExpressionStatement(updateIndexExp);
                //reverse the update
                Expression updateIndexRevExp = new AssignmentExpression(index.clone(), AssignmentOperator.NORMAL, new BinaryExpression(index.clone(), BinaryOperator.DIVIDE, new IntegerLiteral(stepValue)));
                Statement  updateIndexRevExpStat = new ExpressionStatement(updateIndexRevExp);
                Statement loopBody = loop.getBody();
                //add updateIndexExpStat at the begning of the loop
                CompoundStatement newLoopBody = new CompoundStatement();
                newLoopBody.addStatement(updateIndexExpStat.clone());
                newLoopBody.addStatement(loopBody.clone());
                newLoopBody.addStatement(updateIndexRevExpStat.clone());
                loop.setBody(newLoopBody);
                */
                result = true;
            }else{
                result= false;
            }
            
        }
        
       
        
        return result;
    }
    
    public static boolean normalizeLoopInitial(ForLoop loop){
        boolean result = false;
        //get the lrb & upb
        Expression lrb = findLoopLowerBound(loop);
        Expression upb = FindLoopUpperBound(loop);
        if(lrb instanceof IntegerLiteral){
           int lrbValue = (int) ((IntegerLiteral)lrb).getValue();
            //check if the lower bound is 0 
           if(lrbValue == 0){
               //if it is zero no need to transform the loop - just set the result true - which means normalized
               result = true;
           }else if(lrbValue>0){
               //if it is >0 restructure the loop
               //1.set the loop initial statmet as index = 0;
               Expression initExpre = new AssignmentExpression(FindLoopIndex(loop).clone(), AssignmentOperator.NORMAL, new IntegerLiteral(0));
               loop.setInitialStatement(new ExpressionStatement(initExpre.clone()));
               //2. set the uper bound as , upb = upb-lrb
               Expression upbNew = new BinaryExpression(upb.clone(), BinaryOperator.SUBTRACT, lrb.clone());
               //simplify expression
               upbNew = Symbolic.simplify(upbNew);
               loop=setLoopUpperBound(loop, upbNew);
               //3. leave the loop step as itis
               
               //4 replace loop index i by i+lrb  inside the loop body -- reversing the change made at the initial - which is subtracting lrb from index
                Expression oldexpres = FindLoopIndex(loop);
                Expression newExp = new BinaryExpression(oldexpres.clone(), BinaryOperator.ADD, lrb.clone());
                FindAndReplaceExpressions(oldexpres, newExp, loop.getBody());
              
             /*  //old 4 update the loop index at the begining and ending of the loop body
                //index = index+lrb    && index = index-lrb
               Expression indexUpdateExp = new AssignmentExpression(FindLoopIndex(loop).clone(), AssignmentOperator.NORMAL, new BinaryExpression(FindLoopIndex(loop).clone(), BinaryOperator.ADD, lrb.clone()));
               Expression indexUpdateRevExp = new AssignmentExpression(FindLoopIndex(loop).clone(), AssignmentOperator.NORMAL, new BinaryExpression(FindLoopIndex(loop).clone(), BinaryOperator.SUBTRACT, lrb.clone()));
               Statement  indexUpdateStat = new ExpressionStatement(indexUpdateExp.clone());
               Statement  indexUpdateRevStat = new ExpressionStatement(indexUpdateRevExp.clone());
               
               CompoundStatement newLoopBody = new CompoundStatement();
               newLoopBody.addStatement(indexUpdateStat.clone());
               newLoopBody.addStatement(loop.getBody().clone());
               newLoopBody.addStatement(indexUpdateRevStat.clone());
               
               loop.setBody(newLoopBody);
              */ 
               result = true;
           }
           
        }
        return result;
    }
    
    
    /**
     * return the loop step increament value
     * if the step expression is other than increment expression or other than (i=i OP step ) where OP = +; OP is always addition
     * return -1 = which means cant get the step value - and cant be normalized
     * @param loop
     * @return 
     */
    public static int getLoopStepValue(ForLoop loop){
        int value = -1; 
        
        Expression stepExp = loop.getStep();
        if(stepExp instanceof UnaryExpression){
            //get the operator
            UnaryOperator uop =((UnaryExpression)stepExp).getOperator();
            if(uop== UnaryOperator.POST_INCREMENT || uop == UnaryOperator.PRE_INCREMENT){
                value =1;
            }
        }else if(stepExp instanceof AssignmentExpression){
           // System.out.println("Assignment expression");
            //get the right hand expression; since the lesft one is the loop index
             Expression rhs= ((AssignmentExpression)stepExp).getRHS();
            if(rhs instanceof BinaryExpression){
                //check if the operator is ADD
                BinaryOperator bop = ((BinaryExpression)rhs).getOperator();
                if(bop == BinaryOperator.ADD){ //if it is not ADD -- return -1 which means cant find stepvalue
                    //determine if the stepvalue is in the left or right
                   Expression lhs =((BinaryExpression)rhs).getLHS();
                   Expression index = FindLoopIndex(loop);
                   //System.out.println("Add op"+lhs+index);
                   if(lhs.equals(index)){
                      // System.out.println("right side");
                       //then the step value is in the right side
                       Expression rhs2 =((BinaryExpression)rhs).getRHS(); 
                       //check if it is instance of integer literal expression
                       if(rhs2 instanceof IntegerLiteral){
                           value =(int) ((IntegerLiteral)rhs2).getValue();
                       }
                       
                   }else{
                       //the step value is in lhs
                       //check if it is integer literal value
                       if(lhs instanceof IntegerLiteral){
                           value =(int) ((IntegerLiteral)lhs).getValue();
                       }
                   }
                    
                }
            }
        }
        return value;
    }
    
    //FindLoopIndex given a for Loop
    public static Expression FindLoopIndex(ForLoop L){
        Statement stmt = L.getInitialStatement();
        ExpressionStatement es = (ExpressionStatement)stmt;
        Expression expr = es.getExpression();
        BinaryExpression bexpr = (BinaryExpression)expr;
        Expression index = bexpr.getLHS();
        return index;
    }
    
    /**
     * 
     * @param L
     * @return the lower bound of a given forloop
     */
    public static Expression findLoopLowerBound(ForLoop L){
        Statement stmt = L.getInitialStatement();
        ExpressionStatement es = (ExpressionStatement)stmt;
        Expression expr = es.getExpression();
        BinaryExpression bexpr = (BinaryExpression)expr;
        Expression lrb = bexpr.getRHS();
        return lrb;
    }
    
    /**
     * return the forloop upperbound expression
     * assumes that the loop condition expression is binary expression 
     * TODO - check if the loopcondition is inclusive i.e <= if so add 1 to the upperbound - since later the condition is going
     *  to change to only <
     * if not it returns null;
     * @param L
     * @return 
     */
    public static Expression FindLoopUpperBound(ForLoop L){
        Expression upperBound = null;
        Expression condition = L.getCondition();
        if(condition instanceof BinaryExpression){
            BinaryExpression bexpr = (BinaryExpression)condition;
            //if the loop index is at the left- then the upper bound is at the right
            if(bexpr.getLHS().equals(FindLoopIndex(L)))
                upperBound=bexpr.getRHS();
            else
                upperBound= bexpr.getRHS();
        }

        return upperBound;
    }
    
    //find and set the given forloop upperbound expression
    //
    public static ForLoop setLoopUpperBound(ForLoop loop, Expression upb){
        Expression condition = loop.getCondition();
        if(condition instanceof BinaryExpression){
            BinaryExpression bexpr = (BinaryExpression)condition;
            //if the loop index is at the left- then the upper bound is at the right
            if(bexpr.getLHS().equals(FindLoopIndex(loop)))
                bexpr.setRHS(upb);
            else
                 bexpr.setLHS(upb);
        }
        return loop;
    }
    
    //add new function declaration and identifiers
    public static void AddProcedure(Traversable program,Specifier returnSpec,String procName,List params,CompoundStatement body){
        while (!(program instanceof SymbolTable)) {
            program = program.getParent();
        }
        SymbolTable symTable = (SymbolTable)program;
        NameID procNameId= new NameID(procName);
        //create procedure
        Declarator pd = new ProcedureDeclarator(procNameId, params);
        Procedure proc = new Procedure(returnSpec,pd,body);
        
        symTable.addDeclaration(proc);              
        
    }
    
    //Create bool isIterationMine(int rank, int iterationIndex) - function / and add it to program
    //not completed 
    public static void AddComputationPartitioningProcedure(Program program){
              //proc name
           String prcName = "isIterationMine";
                //list of parameters
           List<Declaration> params = new ArrayList<>();
           Declarator par1_dec =new VariableDeclarator(new NameID("rank"));
           Declaration par1 = new VariableDeclaration(Specifier.INT,par1_dec);
           params.add(par1);
           Declarator par2_dec =new VariableDeclarator(new NameID("iterationIndex"));
           Declaration par2 = new VariableDeclaration(Specifier.INT,par2_dec);
           params.add(par2);
                //proc return specifier
           Specifier rtrnSpec = Specifier.BOOL;
                //construct the body 
           CompoundStatement body = new CompoundStatement();
           
           
           MPI_Utility.AddProcedure(MPI_Utility.getMainProcedure(program).getParent(),rtrnSpec, prcName, params,body);
           
    }
    
    //range analysis
    public static void loopRangeAnalysis(ForLoop loop){
        Map<Statement,RangeDomain> range_dom = RangeAnalysis.getRanges(loop,RangeAnalysis.RANGE_INTER);
        Statement body = loop.getBody();
        List<Traversable> stats = body.getChildren();
        for (Traversable stat : stats) {
            RangeDomain rangDo = range_dom.get(stat);
            Expression indexEx = new BinaryExpression(MPI_Utility.FindLoopIndex(loop).clone(), BinaryOperator.SUBTRACT, new IntegerLiteral(3));
            
            //rangDo.replaceSymbol(SymbolTools.getSymbolOf(FindLoopIndex(loop).clone()), indexEx);
           // System.out.println("Statment : "+stat);
            //System.out.println("Range:"+rangDo);
            Expression range =rangDo.getRange(SymbolTools.getSymbolOf(FindLoopIndex(loop)));
            System.out.println("range Exp="+range);
            System.out.println("rangeExpChilds"+range.getChildren().size());
            System.out.println("rangeExpChild0 "+range.getChildren().get(0));
        }
    }
    
    //analyze and parse internal annotation - to find cetus parallel loops
    public static List<ForLoop> getCetusParallelLoops(Program p){
     
        List<ForLoop> parallelLoops = new ArrayList<>();
        //check those inside the main procedure
        //Procedure mainp = getMainProcedure(p);
         /*for each loop in proc find loops*/
        DepthFirstIterator bfs_itr_lp = new DepthFirstIterator(getMainProcedure(p));
           while (bfs_itr_lp.hasNext()) {
               //get the 
               Object o = bfs_itr_lp.next();
               if (o instanceof ForLoop) {
                   ForLoop loop = (ForLoop)o;       
                   
                  //CetusAnnotation cePar = loop.getAnnotation(CetusAnnotation.class, "parallel");
                   //boolean isLoopParallel =IsLoopParallel(dpg,loop);
                   boolean isLoopParallel=true;
                   if(loop.containsAnnotation(CetusAnnotation.class, "parallel") && isLoopParallel){
                       parallelLoops.add(loop);
                       loop.annotate(new CommentAnnotation("Is Parallel="+isLoopParallel));
                   }
                       
                   
                   
               }
           }  
        
       return parallelLoops;
    }
    
    //get reduction variables within the given forloop
    public static void getReductionVariables(ForLoop floop){
        CetusAnnotation ceRedV = floop.getAnnotation(CetusAnnotation.class, "reduction");
        if(ceRedV!=null)
            System.out.println(ceRedV.toString());
    }
    
}
