#include<stdbool.h>

//numProc = size 
//upb = loop upper bound
//lrb = loop loer bound
//rank = processor id
//assuming the loop is normalized, begining index any, stepInc by 1 and loop condition is only < not <=
bool isIterationMine(int numProc,int lrb,int upb, int rank, int iterationIndex){
	bool result = false;
	//given the number of proc, loop upper bound, the processor rank and the iteration index
	//determine whether the given iteration belongs to the given processor
	int blockSize = (upb-lrb)/numProc;
	int blockSizeRemains = (upb-lrb) % numProc;
	int startingIndx; //lower index of the partition for the given processor
	int finishingIndx;  // upper index
	//determine the partition boundary
	if(rank==0){
		startingIndx = lrb;
		finishingIndx = lrb + blockSize+blockSizeRemains;
	}else{
		startingIndx = lrb+(blockSize*rank)+blockSizeRemains;
		finishingIndx = startingIndx + blockSize;
	}
	//check if the given iteration is with in the partition range
	if(iterationIndex >= startingIndx && iterationIndex < finishingIndx)
		result = true;

	return result;

}

//retrun the partition startting point of a processor
//given the number of processors= size; processor rank and acces range [lb:up) -- upper bound is not includeded
//
int getPartitionStartingPoint(int numProc,int lrb,int upb,int rank){
	int partSp;
	int blockSize = (upb-lrb)/numProc;
	int blockSizeRemains = (upb-lrb) % numProc;
	int startingIndx; //lower index of the partition for the given processor
	int finishingIndx;  // upper index
	//determine the partition boundary
	if(rank==0){
		startingIndx = lrb;
		finishingIndx = lrb + blockSize+blockSizeRemains;
	}else{
		startingIndx = lrb+(blockSize*rank)+blockSizeRemains;
		finishingIndx = startingIndx + blockSize;
	}
	partSp=startingIndx;
	return partSp;
}

//return elements of a shared variable that the given node needs
//for now just return the whole element [sv_lb,sv_ub) except those in [p_lb,p_ub) - since they are already owned by the node p
//removing means making the elts -1 ==the elts represent the array subscript

void getElements2send(int sv_lb,int sv_ub,int p_lb, int p_ub, int elements2send[] ){
	int elts2sendSize = sv_ub-sv_lb;
	int eltsNot2sendSize = p_ub-p_lb;
	int i,j;
	//init elements2send - if elts2sendSize > sv_ub-sv_lb -- init the rest with -1 
	//elts2sendSize - the variable size , sv_ub-sv_lb = the variable range that is modified and should be send
	for(i=0;i<elts2sendSize;i++){
		int elt= sv_lb+i;
		//if(elt>=sv_ub) {elt=-1;}
		elements2send[i]=elt;
	}
//remove not2sendelements = [p_lb,p_ub)
	for(i=0;i<elts2sendSize;i++){
		for(j=0;j<eltsNot2sendSize;j++)
		{
		int eltN2send = p_lb+j;
			if(elements2send[i]==eltN2send){
				elements2send[i]=-1;
			}
		}
	}


	
}
