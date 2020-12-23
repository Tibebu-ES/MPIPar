/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mpipar;

import cetus.hir.Expression;
import cetus.hir.PragmaAnnotation;
import cetus.hir.PrintTools;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author tibebu
 */
public class MPI_Annotation extends PragmaAnnotation{
    public static String KEY_ACCESSED_DATA="accessed_data";
    //value indicates loop depth - if 0 = outer loop , if 1=inner loop(inside one loop) 
    public static String KEY_PARALLEL="parallel"; 
    public static String KEY_NORMALIZED="normalized";
    public static String KEY_REDUCTION="reduction";
    public static String KEY_NO_SECOND_COMM="nsc";
    
    public MPI_Annotation(){
        super("MPI");
    }
    
    public MPI_Annotation(String key,Object value){
        super("MPI");
        put(key, value);
    }

    @Override
    public String toString() {
       if (skip_print) {
            return "";
        }
        StringBuilder str = new StringBuilder(80);
        str.append(super.toString()).append(" ");
        if (containsKey(KEY_PARALLEL)) {
            str.append(KEY_PARALLEL).append(" ").append(get(KEY_PARALLEL));
        }
        if (containsKey(KEY_NORMALIZED)) {
            str.append(KEY_NORMALIZED).append(" ");
        }
        if (containsKey(KEY_ACCESSED_DATA)) {
            str.append(KEY_ACCESSED_DATA).append("(");
            List<MPI_AccessedData> sharedVariables =get(KEY_ACCESSED_DATA);
            for (MPI_AccessedData sharedVariable : sharedVariables) {
                if(sharedVariable.isArray())
                    str.append(sharedVariable.getVarExpression()).append(",");
                else
                    str.append(sharedVariable.getVariableSymbol().getSymbolName()+",");
            }
            str.deleteCharAt(str.lastIndexOf(","));
            str.append(")");
        }
        if (containsKey(KEY_REDUCTION)) {
                Map<String, Set<Expression>> reduction_map = this.get(KEY_REDUCTION);
                for (String op : reduction_map.keySet()) {
                    str.append("reduction(").append(op).append(": ");
                    str.append(PrintTools.collectionToString(
                            reduction_map.get(op), ", "));
                    str.append(") ");
                }
        }
        if(containsKey(KEY_NO_SECOND_COMM)){
            str.append(KEY_NO_SECOND_COMM).append(" ");
        }
        return str.toString();
    }
    
    
}
