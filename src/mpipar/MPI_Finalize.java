/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mpipar;

import cetus.hir.Annotation;
import cetus.hir.AnnotationDeclaration;
import cetus.hir.ArraySpecifier;
import cetus.hir.CodeAnnotation;
import cetus.hir.CommentAnnotation;
import cetus.hir.CompoundStatement;
import cetus.hir.Declaration;
import cetus.hir.DeclarationStatement;
import cetus.hir.Declarator;
import cetus.hir.ExpressionStatement;
import cetus.hir.ForLoop;
import cetus.hir.FunctionCall;
import cetus.hir.IDExpression;
import cetus.hir.Identifier;
import cetus.hir.IfStatement;
import cetus.hir.NameID;
import cetus.hir.Procedure;
import cetus.hir.Program;
import cetus.hir.ReturnStatement;
import cetus.hir.Specifier;
import cetus.hir.Statement;
import cetus.hir.SymbolTable;
import cetus.hir.SymbolTools;
import cetus.hir.TranslationUnit;
import cetus.hir.Traversable;
import cetus.hir.UserSpecifier;
import cetus.hir.VariableDeclaration;
import cetus.hir.VariableDeclarator;
import cetus.transforms.LoopNormalization;
import cetus.transforms.TransformPass;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author tibebu
 */
public class MPI_Finalize  extends TransformPass{

    public MPI_Finalize(Program program){
        super(program);
       
    }
    
    @Override
    public String getPassName() {
        return "MPI Finalize.";
    }

    @Override
    public void start() {
        //System.out.println("disable protection ="+disable_protection+"disable invalidation "+disable_invalidation);
        addHeaderFiles(program);
        transformMainProcParams(program);
        addMPIVariablesFunctions(program);
        
       // setOutputFileName(program); 
           
    }
    
    /**
     * add the include header files at the top of source file (mpi.h mpi_utility.h and string.h)
     * and create mpi_utility.h translation unit and add it to the program 
     * @param p 
     */
    public  void addHeaderFiles(Program p){
        //assume the program has one translation unit
        //get the program trnslation unit: program  may contain many source files - where each are representedby translation unit
        TranslationUnit tu = MPI_Utility.getTranslationUnit(p);
        
        //create the annotation declaration
        String rawCode = "#include \"mpi.h\" \n"
                + "#include <string.h> \n"
                + "#define TAG 123";
        Annotation anno1 = new CodeAnnotation(rawCode);
        AnnotationDeclaration anno1Dec = new AnnotationDeclaration(anno1);
        
        //add declaration at the head of the translation unit
        tu.addDeclarationFirst(anno1Dec);
     /*   
        //create & add the "utility.h" source file - as a translation unit to the program
        //1st create utility.h translation unit - from the utility.h source file in the myresource folder
        TranslationUnit utilityTU = sourceFileToTranslationUnit("mpi_utility.h");
        //then add this translation unit to the program
        if(utilityTU!=null)
            p.addTranslationUnit(utilityTU);
     */   
        
    }
    
    /**
     * Given the name of the source file in the myheaders folder - 
     * it returns a translation unit instance - out of the given source file
     * returns null if the file can be find or read 
     * @param fileName
     * @return 
     */
    public  TranslationUnit sourceFileToTranslationUnit(String fileName){
        TranslationUnit tu = new TranslationUnit(fileName);
        try{
  
             //DataInputStream di = new DataInputStream(getClass().getClassLoader().getResourceAsStream("myheaders/"+fileName));
             InputStreamReader isr = new InputStreamReader(getClass().getClassLoader().getResourceAsStream("myheaders/"+fileName));
             BufferedReader brd = new BufferedReader(isr);
             String sourceCode ="";
             int i;
             while((i= brd.read())!=-1)
                 sourceCode+=(char)i;
             
             brd.close();
             
             //Now create annotation declaration from the source code
             AnnotationDeclaration ad = new AnnotationDeclaration(new CodeAnnotation(sourceCode));
             //add this declaration to the translation unit
             tu.addDeclaration(ad);
         
        }catch(Exception e){
           tu = null;
            System.out.println("Error:"+e.toString());
        }
        
        return tu;
    }
    /**
     * check if the required parameters of the main procedure already exist
     * required parameters are (int argc, char *argv[])
     * if not exist - add the parameters 
     * @param p 
     */
    public  void transformMainProcParams(Program p){
        //get the main procedure
        Procedure main = MPI_Utility.getMainProcedure(p);
        
        
        //declare the parameters  - i.e  int argc , char *argv[]
        Declarator par1_dec =new VariableDeclarator(new NameID("argc"));
        Declaration par1 = new VariableDeclaration(Specifier.INT,par1_dec);
        //ArraySpecifier
        Declarator par2_dec =new VariableDeclarator(new NameID("*argv"),ArraySpecifier.UNBOUNDED);
        List<Specifier> par2_specs = new ArrayList<>();
        par2_specs.add(Specifier.CHAR);
        //par2_specs.add(ArraySpecifier.UNBOUNDED);
        Declaration par2 = new VariableDeclaration(par2_specs,par2_dec);
        //check if they exist already
        List<VariableDeclaration> params = main.getParameters();
        boolean containsIntArgcParam = false;
        boolean containsCharArgvParam = false;
        for (VariableDeclaration param : params) {
            //if(param.getDeclarator(0).getID().toString().equals(par1_dec.getID().toString()))
            IDExpression varName = param.getDeclaredIDs().get(0);
            List<Specifier> varSpecifiers = param.getSpecifiers();
            if(varName.equals(par1_dec.getID()) && varSpecifiers.contains(Specifier.INT))
                containsIntArgcParam=true;
            //Todo- check for the [] & * specifiers|
            if(varName.equals(new NameID("argv")) && varSpecifiers.contains(Specifier.CHAR))
                containsCharArgvParam=true;
            
            if(containsCharArgvParam && containsIntArgcParam )
                break;
                
        }
        if(!containsIntArgcParam)
            main.addDeclaration(par1);
        if(!containsCharArgvParam)
            main.addDeclaration(par2);
    }
    
    
    /**
     * Insert necessary Mpi variables with the program main method like
     * MPI_Status status, int rank , int size.
     * add basic mpi functions i.e
     * MPI_Init(&argc,&argv);MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD, &size);
     * and Add MPI_Finalize(); function before the return statment of the main procedure.
     * @param p 
     */
    public  void addMPIVariablesFunctions(Program p){
        Traversable mainBody = MPI_Utility.getMainProcedure(p).getBody();
        //find scope
        while (!(mainBody instanceof SymbolTable)) {
            mainBody = mainBody.getParent();
        }
        SymbolTable symtab = (SymbolTable) mainBody;
        
        //add comment annotation - MPI Variables declaration and MPI basic Function call begin
        Annotation mpiVariablesBeginComment = new CommentAnnotation("MPI variables declaration  begins");
        Annotation mpiVariablesEndComment = new CommentAnnotation("MPI variables declaration ends");
        
        //MPI variable  - MPI_Status status
        Declarator mpiStatusD = new VariableDeclarator(new NameID("status"));
        Declaration mpiStatusDec = new VariableDeclaration(new UserSpecifier(new NameID("MPI_Status")),mpiStatusD);
        //add comment annotation - MPI Variables declaration begin
        mpiStatusDec.annotateBefore(mpiVariablesBeginComment);
        symtab.addDeclaration(mpiStatusDec);
        
        //MPI variable  - int rank
        Declarator rankD = new VariableDeclarator(new NameID("rank"));
        Declaration rankDec = new VariableDeclaration(Specifier.INT,rankD);
        symtab.addDeclaration(rankDec);
        
        //Declare other mpi variables here
        
        
        
         //MPI variable  - int size
        Declarator sizeD = new VariableDeclarator(new NameID("size"));
        Declaration sizeDec = new VariableDeclaration(Specifier.INT,sizeD);
       // sizeDec.annotateAfter(mpiVariablesEndComment);
        symtab.addDeclaration(sizeDec);
        
        //Add MPI basic function - MPI_Init(&argc,&argv);
        Identifier mpi_init = SymbolTools.getOrphanID("MPI_Init");
        FunctionCall mpi_initFcnCall = new FunctionCall(mpi_init.clone());
        mpi_initFcnCall.addArgument(new NameID("&argc"));
        mpi_initFcnCall.addArgument(new NameID("&argv"));
        Statement mpi_initFcnCallStmnt = new ExpressionStatement(mpi_initFcnCall);
        
        //add it after the last mpi variable declaration statment
        //get the reference statment which is the last mpi variable declaration
        Statement lastMPIVariableDecStat = null;
        CompoundStatement mainBC = MPI_Utility.getMainProcedure(p).getBody();
        List<Traversable> childerns = mainBody.getChildren();
        for (Traversable childern : childerns) {
            if(childern instanceof DeclarationStatement){
                DeclarationStatement s = (DeclarationStatement)childern;
                if(s.getDeclaration().equals(sizeDec)){
                    lastMPIVariableDecStat=s;
                    break;
                }
            }
        }
        
        mainBC.addStatementAfter(lastMPIVariableDecStat, mpi_initFcnCallStmnt);
        
        //Add MPI basic function - MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        Identifier MPI_Comm_rank = SymbolTools.getOrphanID("MPI_Comm_rank");
        FunctionCall MPI_Comm_rankFcnCall = new FunctionCall(MPI_Comm_rank.clone());
        MPI_Comm_rankFcnCall.addArgument(new NameID("MPI_COMM_WORLD"));
        MPI_Comm_rankFcnCall.addArgument(new NameID("&rank"));
        Statement MPI_Comm_rankFcnCallStmnt = new ExpressionStatement(MPI_Comm_rankFcnCall);
        mainBC.addStatementAfter(mpi_initFcnCallStmnt,MPI_Comm_rankFcnCallStmnt);
        
        //Add MPI basic function - MPI_Comm_size(MPI_COMM_WORLD, &size);
        Identifier MPI_Comm_size = SymbolTools.getOrphanID("MPI_Comm_size");
        FunctionCall MPI_Comm_sizeFcnCall = new FunctionCall(MPI_Comm_size.clone());
        MPI_Comm_sizeFcnCall.addArgument(new NameID("MPI_COMM_WORLD"));
        MPI_Comm_sizeFcnCall.addArgument(new NameID("&size"));
        Statement MPI_Comm_sizeFcnCallStmnt = new ExpressionStatement(MPI_Comm_sizeFcnCall);
        //add comment annotation - MPI Variables declaration begin
        MPI_Comm_sizeFcnCallStmnt.annotateAfter(mpiVariablesEndComment);
        mainBC.addStatementAfter(MPI_Comm_rankFcnCallStmnt,MPI_Comm_sizeFcnCallStmnt);
        
       
        //Add MPI_Finalize(); function before the return statment of the main procedure.
        FunctionCall MPI_Finalize = new FunctionCall(SymbolTools.getOrphanID("MPI_Finalize"));
        Statement MPI_FinalizeStmnt = new ExpressionStatement(MPI_Finalize);
        //get the return statement
        Statement retrurnStmnt = null;
        for (Traversable childern : childerns) {
            if(childern instanceof ReturnStatement){
                retrurnStmnt = (Statement) childern;
                break;
            }
        }
        
        if(retrurnStmnt!=null){
            mainBC.addStatementBefore(retrurnStmnt, MPI_FinalizeStmnt);
        }
        
    }
    

    /**
     * set the output file name to inputName_mpi.c
     * simply add '_mpi.c' to the name of the input file
     */
    public void setOutputFileName(Program program){
        //change the the name of output file
        TranslationUnit mainTu = MPI_Utility.getTranslationUnit(program);
        String originalOutFname = mainTu.getOutputFilename();
        // System.out.println("TU name="+originalOutFname);
        mainTu.setOutputFilename(originalOutFname.substring(0, originalOutFname.indexOf('.')) + "_mpi.c");
        //System.out.println("TU name="+mainTu.getOutputFilename());
    }

}

