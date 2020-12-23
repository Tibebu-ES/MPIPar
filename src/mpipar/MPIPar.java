/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mpipar;

import cetus.analysis.AliasAnalysis;
import cetus.analysis.AnalysisPass;
import cetus.analysis.Domain;
import cetus.analysis.IPAnalysis;
import cetus.analysis.IPPointsToAnalysis;
import cetus.analysis.PointsToAnalysis;
import cetus.exec.Driver;
import static cetus.exec.Driver.getOptionValue;
import cetus.hir.BreadthFirstIterator;
import cetus.hir.Declaration;
import cetus.hir.IDExpression;
import cetus.hir.Identifier;
import cetus.hir.NameID;
import cetus.hir.Statement;
import cetus.hir.Symbol;
import cetus.hir.SymbolTools;
import cetus.hir.Traversable;
import cetus.hir.VariableDeclaration;
import cetus.transforms.TransformPass;
import java.util.Set;

/**
 *
 * @author tibebu
 */
public class MPIPar extends Driver {

    public String comm_scheme_selected;
   
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        // TODO code application logic here
        if (args == null || args.length == 0) {
            //"-tinline=mode=2:functions=kernel_2mm"
            String myargs[] = {"-range=1","-alias=1","-outdir=mpi_output", "-MPIGen", "-TIME","mvt.c"};
            args = myargs;
        }

        (new MPIPar()).run(args);
    }

    public MPIPar() {
        super();
       
        options.add("MPIGen", "Generate MPI code");
        options.add("TIME", "Insert comm overhead instrumenting codes");
    }

    public String getCommSchemeSelected(){
        return  this.comm_scheme_selected;
    }
    
    @Override
    public void runPasses() {
        super.runPasses(); //To change body of generated methods, choose Tools | Templates.

        //reading command line argument
        boolean doMPI = false;
        boolean timeReport = false;
        String valueMPIGen = getOptionValue("MPIGen");
        String valueTIME = getOptionValue("TIME");
        if (valueMPIGen != null) {
            doMPI = true;
        }
        if(valueTIME!=null)
            timeReport=true;

        
        //run my pass
        if (doMPI) {
         //try to parse the use choice of comm schem and set the value
             
           // try{
                 //System.out.println("My  MPIGen to goes here!");
                //run MPI normalize trnasform
                TransformPass.run(new MPI_LoopNormalize(program));
           //run MPI analysis - find cetus paralle loops - mpi normalized loops and attach mpi parallel annotation
                //also find shared variables in each parallel loop - attach MPI sharedVariables(list of variables) annotation
                AnalysisPass.run(new MPI_Analysis(program));
                //run MPI transform - add mpi constructs
               TransformPass.run(new MPI_Finalize(program));

                //--detach annotations except the MPI_Annotations parallel & sharedvariables- loop partition - communication code generation
                //set the default comm scheme
               //use default comm scheme
              TransformPass.run(new MPI_LoopTransform(program,MPI_LoopTransform.COMM_SCHEME_ONE,timeReport));
         //   }catch(Exception e){
         //       System.out.println("MPIGen Error :"+e.getMessage());
         //   }
            

        }
    }

}
