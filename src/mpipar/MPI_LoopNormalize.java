/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mpipar;

import cetus.analysis.Domain;
import cetus.analysis.IPPointsToAnalysis;
import cetus.analysis.LoopTools;
import cetus.hir.Annotation;
import cetus.hir.CetusAnnotation;
import cetus.hir.CommentAnnotation;
import cetus.hir.ForLoop;
import cetus.hir.Loop;
import cetus.hir.PragmaAnnotation;
import cetus.hir.Program;
import cetus.transforms.LoopTransformPass;
import java.util.List;

/**
 *try to normalize all for loops 
 * add MPI normalized -- annotation to those normalized correctly
 * 
 * which will be later used by the MPI_Analysis pass 
 * @author tibebu
 */
public class MPI_LoopNormalize extends LoopTransformPass{
    
    public MPI_LoopNormalize(Program program){
        super(program);
         //disable analysis consistency checking
          disable_protection = true;
          disable_invalidation=true;
    }

    @Override
    public void transformLoop(Loop arg0) {
        if(arg0 instanceof ForLoop){
            ForLoop floop = (ForLoop) arg0;
        //try ton normalize
            boolean isNormalized = MPI_Utility.NormalizeLoop(floop);
            //boolean isNormalized = true;
            if(isNormalized){
                //attach MPI normalized annotation
                MPI_Annotation mpiNormalizedAn = new MPI_Annotation("normalized", 1);
                floop.annotate(mpiNormalizedAn);
            }else
                floop.annotate(new CommentAnnotation("Not Normalized"));
       /*     
            //remove all annotations of loops outside the main pro
            if(!floop.getProcedure().equals(MPI_Utility.getMainProcedure(program))){
                floop.removeAnnotations();
            }else{
                
                //for loops in the main proc
                //deattach annotations except cetusAnnotation & MPI annotations & 
                List<CetusAnnotation> cetAnno =floop.getAnnotations(CetusAnnotation.class);
                List<MPI_Annotation> mpiAnno =floop.getAnnotations(MPI_Annotation.class);
                floop.removeAnnotations();

                if(!cetAnno.isEmpty()){
                    for (CetusAnnotation cetAnno1 : cetAnno) {
                        floop.annotate(cetAnno1);
                    }
                }
                if(!mpiAnno.isEmpty()){
                    for (MPI_Annotation mpiAnno1 : mpiAnno) {
                        floop.annotate(mpiAnno1);
                    }
                }
                
            }
         */   
           //deattaching annotations finished    
        }
    }

    @Override
    public String getPassName() {
        return "MPI Loop Normalize";
    }
    
    public void aliasAnalysis(ForLoop floop){
                if(LoopTools.getLoopName(floop).contains("main")){
                 Domain dm=IPPointsToAnalysis.getPointsToRelations(floop);
                 System.out.println(LoopTools.getLoopName(floop)); 
                 System.out.println(dm.toString()); 
                 //iterate through each statments
                 
                 
            }
         
        
    }
    
}
