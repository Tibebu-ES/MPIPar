/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mpipar;

import cetus.analysis.RangeAnalysis;
import cetus.analysis.RangeDomain;
import cetus.hir.ArrayAccess;
import cetus.hir.ArraySpecifier;
import cetus.hir.BinaryExpression;
import cetus.hir.BinaryOperator;
import cetus.hir.Expression;
import cetus.hir.IntegerLiteral;
import cetus.hir.NameID;
import cetus.hir.PointerSpecifier;
import cetus.hir.RangeExpression;
import cetus.hir.Specifier;
import cetus.hir.Symbol;
import cetus.hir.SymbolTable;
import cetus.hir.SymbolTools;
import cetus.hir.Symbolic;
import cetus.hir.Traversable;
import cetus.hir.TypeofSpecifier;
import java.util.ArrayList;
import java.util.List;


/**
 * contains information about the accessed data neccessary to generate communication code like
 * if it is scalar or array
 * the variable symbol
 * the datatype
 * if reduction variable or not 
 * index of the parallel loop, where this accessed data belongs
 *  for Array type
 *   - array access
 *   - number of dimensions
 *   - dimension where the parallel loop index is found
 *   - block multiplier - which is number of elements to the right of the dimension where parallel loop is found
 * @author tibebu
 */
public class MPI_AccessedData {
    
    //communication types
    public static String   COMMUNICATION_TYPE_BROADCASTING="broadcasting";
    public static String   COMMUNICATION_TYPE_GATHERING="gathering";
    /**
     * variable symbol
     */
   private Symbol variableSymbol;
    /**
     * expression that contains the accessed variable
     * used if the accessed data is array type
     */
   private Expression varExpression;
   /**
    * represents the access of this accessed data -if it is array - array[x][y]
    */
   private ArrayAccess arrayAccess;
   /**
    * RangeDomain - holds the range of all symbols participated in the index expression of this variable
    * used to calculate the access range of each dimension
    */
   private RangeDomain rangeDomain;
   /**
   * holds info - if shared variable is reduction variable or not
   */
   private boolean isReductionVariable;
   /**
    * reduction operator- for reduction variables
    */
   private BinaryOperator reductionOperator;
   /**
    * holds the index symbol of the parallel loop - where this var belongs
    * */
   private Symbol indexOftheParallelLoop;
   /**
    * holds if it is scalar or array variable
    * */
   private boolean isArray;
   /**
    * holds - type of communication needed - broadcasting or gathering 
    * broadcasting is the default
    */
   private String commType;
   
   
    public MPI_AccessedData(Symbol variableSymbol) {
        this.variableSymbol = variableSymbol;
        
        
        this.reductionOperator = null;
        this.isReductionVariable = false;
        //default comm type is broadcasting
        this.commType=COMMUNICATION_TYPE_BROADCASTING;
    }
    
    
    //setters and getters

    
    public void setCommType(String commType) {
        this.commType = commType;
    }

    public String getCommType() {
        return commType;
    }
    
    

    public void setIsArray(boolean isArray) {
        this.isArray = isArray;
    }

    public boolean isArray() {
        return isArray;
    }
    
    

    public void setIndexOftheParallelLoop(Symbol indexOftheParallelLoop) {
        this.indexOftheParallelLoop = indexOftheParallelLoop;
    }

    public Symbol getIndexOftheParallelLoop() {
        return indexOftheParallelLoop;
    }

    /**
     * if the accessed data is array
     * returns the dimension where  the parallel loop index is found
     * @return 
     */
    public int getDimOfIndexOftheParallelLoop(){
        int dimNum =-1;
        if(getIndexOftheParallelLoop()!=null)
            dimNum=getDimNumOfIndex(getIndexOftheParallelLoop());
        
        return dimNum;
    }

    /**
     * also set isReduction boolean 
     * @param reductionOperator 
     */
    public void setReductionOperator(BinaryOperator reductionOperator) {
        if(reductionOperator!=null){
            this.reductionOperator = reductionOperator;
            this.isReductionVariable=true; 
        }
    }

    public BinaryOperator getReductionOperator() {
        return reductionOperator;
    }

    public boolean isIsReductionVariable() {
        return isReductionVariable;
    }
    
    public void setVariableSymbol(Symbol variableSymbol) {
        this.variableSymbol = variableSymbol;
    }

    public Symbol getVariableSymbol() {
        return variableSymbol;
    }

    public Specifier getDataType() {
        Specifier dataType =null;
        if(getVariableSymbol()!=null){
            dataType = (Specifier)getVariableSymbol().getTypeSpecifiers().get(0);
        }
        return dataType;
    }

    
    /**
     * sets the varExpression which- triggers the varExpressionToArrayAccess and set the array access attribute
     * @param varExpression 
     */
    public void setVarExpression(Expression varExpression) {
        this.varExpression = varExpression;
        setArrayAccess(arrayExpressionToArrayAccess(varExpression));
    }

    public Expression getVarExpression() {
        return varExpression;
    }
    
    /**
     * 
     * @return the dimension of this variableArray 
     */
    public int getNumOfDimension(){
        ArraySpecifier as =(ArraySpecifier) this.variableSymbol.getArraySpecifiers().get(0);
      /*  int numOfDim = 0;
        if(getArrayAccess()!=null)
            numOfDim=getArrayAccess().getNumIndices();
        return numOfDim;
        */
        return as.getNumDimensions();
        
    }
    
    //ABOUT THE ARRAY ACCESS - getting the access range of each dimension
    /**
     * return array access expression of accessed data - if it is array type
     * @return 
     */
    public ArrayAccess getArrayAccess() {
        return arrayAccess;
    }

    public void setArrayAccess(ArrayAccess arrayAccess) {
        this.arrayAccess = arrayAccess;
    }
    
    /**
     * return array access from the given var expression - like array[x][y]
     * @param varExpression
     * @return 
     */
    public ArrayAccess arrayExpressionToArrayAccess(Expression varExpression){
        ArrayAccess arrayAccess = null;
        //iterate through the varexpression children and get the list of each dimension index expression -- starting from index 0
        List<Expression> indices = new ArrayList();  
        List<Traversable> children = varExpression.getChildren();
        //the first child is the array expression
        Expression arrayExp = ((Expression)children.get(0)).clone();
        //start from the second child since the first child is the array name ---the   
        for(int i=1;i<children.size();i++){
               //child other than the first one is expression of the subscripts
            indices.add(((Expression)children.get(i)).clone());
        }
        
        arrayAccess = new ArrayAccess(arrayExp, indices);
        return arrayAccess;
    }

    public void setRangeDomain(RangeDomain rangeDomain) {
        this.rangeDomain = rangeDomain;
        //check substitue forward
        if(this.rangeDomain!=null)
            this.rangeDomain.substituteForwardRange();
    }

    public RangeDomain getRangeDomain() {
        return rangeDomain;
    }
    
    /**
     * calculate access range of the nth dimension
     * if there is problem calculating the range access --- it will return -1:-1 range 
     * returns the acces range as [lrb:up) - which means the upper bound is not included
     * @param n
     * @return rangeExpression(lb,ub)
     * 
     */
    public RangeExpression getAccessRange(int n){
        RangeExpression rangeExpression = new RangeExpression(new IntegerLiteral(-1), new IntegerLiteral(-1));
        //if rangeDom is null or n is greather than the array dimension or  return -1:1
        if(getRangeDomain() ==null || n>=getNumOfDimension() || n<0)
            return rangeExpression;
        //1. get the symbols involved in the index expression
       
        Expression indexExpr = getArrayAccess().getIndex(n);
       // System.out.println("array access"+getArrayAccess().getIndex(n));
        //put it in a list
        List<Expression> indexExprList = new ArrayList<>();
        indexExprList.add(indexExpr);
        List<Symbol> symbols = SymbolTools.exprsToSymbols(indexExprList);
        //2 calculate the access range lower bound and upper bound
        Expression lb = indexExpr.clone();
        Expression ub = indexExpr.clone();
        //System.out.println("indexExpr"+indexExpr);
        //2.1 iterate through childs of expression & check if the child is a symbol - if so
        //find the lower bound and upper bound of that symbol and replace these values by the childs of the lb & up express
        //respectively
        List<Traversable> children =indexExpr.getChildren();
       // System.out.println("expr children num="+children.size());
        //for empty children
        if(children.isEmpty()){
           Symbol indexExpChildSym= SymbolTools.getSymbolOf(indexExpr);
           if(indexExpChildSym != null){
                //if it is symbol - then get the range expression - which is like [lb:ub:stride]
                 //Expression symbolRange=RangeAnalysis.getRangeDomain(indexExpr.getStatement()).getRange(indexExpChildSym);
                //System.out.println("Rangedomain="+getRangeDomain());
                Expression symbolRange = getRangeDomain().getRange(indexExpChildSym);
                System.out.println("SymbolRange of"+indexExpChildSym+"="+symbolRange);
                if(symbolRange!=null){
                     Expression symbol_lb = (Expression) symbolRange.getChildren().get(0); //first child is lb
                    Expression symbol_ub = (Expression) symbolRange.getChildren().get(1); //second child is ub
                    //now replace the corresponding childs of the cloned lb ub expressions
                    lb=symbol_lb.clone();
                    
                    ub=symbol_ub.clone(); 
                }
              
            }
        }
        for(int i=0;i<children.size();i++){
            //check if the child is symbol - if not return value will be null    
            Symbol   indexExpChildSym= SymbolTools.getSymbolOf((Expression) children.get(i));
            //System.out.println("indexExpChildSym"+indexExpChildSym);
            if(indexExpChildSym != null){
                //if it is symbol - then get the range expression - which is like [lb:ub:stride]
                 //Expression symbolRange=RangeAnalysis.getRangeDomain(indexExpr.getStatement()).getRange(indexExpChildSym);
                Expression symbolRange = getRangeDomain().getRange(indexExpChildSym);
                //System.out.println("SymbolRange of"+indexExpChildSym+"="+symbolRange);
                if(symbolRange!=null){
                    Expression symbol_lb = (Expression) symbolRange.getChildren().get(0); //first child is lb
                    Expression symbol_ub = (Expression) symbolRange.getChildren().get(1); //second child is ub
                    //now replace the corresponding childs of the cloned lb ub expressions
                    lb.setChild(i, symbol_lb.clone());
                    ub.setChild(i, symbol_ub.clone()); 
                }
                
            }
        }
        //since the upper shouldnt be included - add one to it
        ub = new BinaryExpression(ub, BinaryOperator.ADD, new IntegerLiteral(1));
        //3 simplify the lb and ub expressions
        lb = Symbolic.simplify(lb);
        ub = Symbolic.simplify(ub);
        
        
        rangeExpression = new RangeExpression(lb, ub);
        return rangeExpression;
    }
    
    /**
     * iterate through the arrayacces indices - find in which dim is the given index located
     * if not found return -1
     * if arrayaceess is null also return -1
     * @param index
     * @return 
     */
    public int getDimNumOfIndex(Symbol index){
        int dimNum = -1; //if it is not found -1 will be returned
        ArrayAccess aa = this.getArrayAccess();
        if(aa!=null){
            for(int dim =0;dim<aa.getNumIndices();dim++){
                Expression indexExp = new NameID(index.getSymbolName());
                if(!aa.getIndex(dim).findExpression(indexExp).isEmpty()){
                    //index find here
                    dimNum=dim;
                    break;
                }
            }
        }
        return dimNum;
    }
    
    /**
     * return the product of size of dimensions below the partitioned dimension
     * returns the number of elements of sv - starting from the
 getDimOfIndexOftheParallelLoop()+1
     * @return 
     */
    public Expression getPartitionSizeMultiplier(){
        Expression mult = new IntegerLiteral(1);
        if(getDimOfIndexOftheParallelLoop()!=-1){
        //System.err.println("From BSM---"+getVariableSymbol().getSymbolName()+"=NUmofDIm="+getNumOfDimension());    
            for(int i=getDimOfIndexOftheParallelLoop()+1;i<getNumOfDimension();i++){
                ArraySpecifier as =(ArraySpecifier) getVariableSymbol().getArraySpecifiers().get(0);
               // System.err.println("From BSM---Decla="+getVariableSymbol().getDeclaration());
                mult = new BinaryExpression(mult.clone(), BinaryOperator.MULTIPLY,as.getDimension(i).clone());
               // System.err.println("From BSM---"+getVariableSymbol().getSymbolName()+"="+mult);
            }
        }
        return Symbolic.simplify(mult);
        //return mult;
    }
    
    /**
     * return the var expression with
     * parallel loop index expression in the array subscript replaced by the given expression &
 dimensions after the getDimOfIndexOftheParallelLoop() - replaced by its lb expression
 use the arrayacces expression 
     * @param partSp
     * @return 
     */
    public Expression getVarExprWithParLpIndexReplaced(Expression partSp){
        ArrayAccess accesExpr = getArrayAccess().clone();
        //Expression partSp = new NameID("partSp");
        if(accesExpr!=null){
            //replace subscript expr of dimensions after the getDimOfIndexOftheParallelLoop() - by its lb expression
            for(int i=getDimOfIndexOftheParallelLoop()+1;i<getNumOfDimension();i++){
                Expression lb = getAccessRange(i).getLB().clone();
                accesExpr.setIndex(i, lb.clone());
            }
            //get the index expression where - the parallel loop index located
           // Expression indexExp =accesExpr.getIndex(getDimOfIndexOftheParallelLoop());
            //now find and replace the parrel loop index by - partSp
            //Utility.FindAndReplaceExpressions(getIndexOftheParallelParent(), partSp, indexExp);
            //replace subscript expr of dimension where indexOftheParallel loop located by partSp
             accesExpr.setIndex(getDimOfIndexOftheParallelLoop(), partSp.clone());
            
        }
     
        return accesExpr;
    }
    
    @Override
    public String toString() {
        StringBuilder str = new StringBuilder();
        str.append("Variable name =").append(this.variableSymbol.getSymbolName()).append("\n");
        str.append("Data type =").append(getDataType().toString()+"\n");
        str.append("isReductionVariable = "+isReductionVariable+"\n");
        str.append("Comm Type ="+getCommType()+"\n");
        if(isArray){
            str.append("Type = Array\n");
            str.append("Array Access ="+this.getArrayAccess()+"\n");
            str.append("Number of Dim ="+getNumOfDimension()+"\n");
            str.append("Range Domain ="+getRangeDomain()+"\n");
            str.append("DimNumofParentLoop("+getIndexOftheParallelLoop()+") ="+getDimOfIndexOftheParallelLoop()+"\n");
            str.append("Access range of dim("+getDimOfIndexOftheParallelLoop()+")="+getAccessRange(getDimOfIndexOftheParallelLoop())+"\n");
            str.append("BlockSize Multiplier ="+getPartitionSizeMultiplier()+"\n");
         
            //str.append("var exp with ="+getVarExprWithParLpIndexReplaced(new NameID("partSp"))+"\n");
            
        }else{
            str.append("Type = Scalar\n");
        }
               
        
        return str.toString();
    }

    @Override
    public boolean equals(Object o) {
        
        boolean isEqual = false;
        if(o instanceof MPI_AccessedData){
            MPI_AccessedData o2 = (MPI_AccessedData)o;
            //for array types
            if(isArray){
                 if(o2.isArray && (getVarExpression().toString().equals(o2.getVarExpression().toString())))
                     isEqual=true;
                 
            }
            else{
                if(!o2.isArray && (getVariableSymbol().getSymbolName().equals(o2.getVariableSymbol().getSymbolName())))
                     isEqual=true;
                
            }
           
        }
       
        
        return isEqual;
       // return super.equals(o);
    }
    
}
