#include "ODFCost.h"

static unsigned int iterationNum = 0; 

// initialize internal variables
ODFCost::ODFCost ( ODFCost::ODFImageType::Pointer odf, DoubleImageType::Pointer origfa, ODFCost::LabelImageType::Pointer wmMask, unsigned long verts ,double alpha) // WM mask added by Adrien Kaiser
{
   this->m_ODFImage = odf ;

   //read FA file
   this->m_OrigFAImage = origfa;

   // WM mask added by Adrien Kaiser
   this->m_WMmaskArray = wmMask->GetPixelContainer()->GetBufferPointer();

   this->m_alpha = alpha;
   this->m_NumberOfDirs = verts ;
   this->m_NumberOfSourceVoxels = 0 ;
   this->m_ODFArray = odf->GetPixelContainer()->GetBufferPointer();

   
   this->m_OrigFAArray = origfa->GetPixelContainer()->GetBufferPointer();

   this->PrecomputeStepLength () ;
   this->FstarTraversalUpdates () ;
   std::cout << "initialization done" << std::endl ;
}      

// main function - public
ODFCost::DoubleImageType::Pointer ODFCost::Cost ( ODFCost::LabelImageType::Pointer source ) 
{
  this->m_Converged = false ;
  this->m_NChanged = 0 ;
  this->m_SourceArray = source->GetPixelContainer()->GetBufferPointer();

  this->CreateCostMaps () ;
  
  this->FindSourceVoxels () ;
  
  this->InitializeSourceVoxels () ;
 
  this->NormalizeODF () ;
  
  this->CreateFinslerTable();

  this->PrecomputeDirection();

  this->PrecomputeNeighbors () ;

//   this->m_NeighborArray = this->m_FinalCostImage->GetPixelContainer()->GetBufferPointer() ;

  std::cout << "Entering main loop" << std::endl ;
  unsigned i = 0 ;

/* Display constant values*/
unsigned long int x1=46;
unsigned long int y1=47;
unsigned long int z1=38;
unsigned long int x2=46;
unsigned long int y2=48;
unsigned long int z2=38;

unsigned long voxelIndex1 = x1 + y1 * this->m_XDim + z1 * this->m_Slice ;
unsigned long voxelIndexTimesNDirs1 = voxelIndex1 * this->m_NumberOfDirs ;
unsigned long voxelIndex2 = x2 + y2 * this->m_XDim + z2 * this->m_Slice ;
unsigned long voxelIndexTimesNDirs2 = voxelIndex2 * this->m_NumberOfDirs ;

for ( unsigned long dir = 0 ; dir < this->m_NumberOfDirs ; dir++ )
{
	unsigned int voxelAndIn1 = voxelIndexTimesNDirs1 + dir ;
	unsigned int voxelAndIn2 = voxelIndexTimesNDirs2 + dir ;
	std::cout<<"dir=,"<<dir<<" \t,"<<this->m_CoordinateTable[dir][0]<<"\t,"<<this->m_CoordinateTable[dir][1]<<"\t,"<<this->m_CoordinateTable[dir][2]<<"\t,]"<<std::endl;
//	std::cout<<"dir="<<dir<<" \t: Finsler1=,"<<this->m_FinslerTable[voxelAndIn1]<<" , Finsler2=,"<<this->m_FinslerTable[voxelAndIn2]<<std::endl;
}

  while ( !this->m_Converged ) 
    {
      iterationNum ++ ;
      this->MakeOnePassFstar (i) ;
      i++ ;
      std::cout << "Pass " << i << " complete. " << this->m_NChanged << " voxels updated in this pass." << std::endl ;
      
      if(this->m_NChanged == 0)
	{
	  this->m_Converged = true;
	}
      else
	{
	  this->m_Converged = false;
	}    
      //---------------------------------------
      this->m_NChanged = 0 ;
    }

  std::cout << "Compute final costs" << std::endl ;
  this->ComputeFinalCosts () ;

  return this->m_FinalCostImage ;
}

// visit every voxel in the image
void ODFCost::MakeOnePass ()
{
  
  unsigned long int x, y, z ;

  for ( z = 0 ; z < this->m_ZDim ; z++ )
    {
      for ( y = 0 ; y < this->m_YDim ; y++ )
	{
	  for ( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      this->VisitVoxel ( x, y, z ,0,0) ;
	    }
	}
   }
}


// visit every voxel according to F* traversal algorithm
void ODFCost::MakeOnePassFstar (int i)
{
  short x, y, z ;

  for( z = 0 ; z < this->m_ZDim ; z++ )
    {
      for( y = 0 ; y < this->m_YDim ; y++ )
	{
	  this->m_fstar_direction = 0 ;
	  for( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      this->VisitVoxel ( x, y, z ,1,i) ;
	    }
	  this->m_fstar_direction = 1 ;
	  for( x = this->m_XDim - 1 ; x >= 0 ; x-- )
	    {
	      this->VisitVoxel ( x, y, z ,2,i) ;
	    }
	}
      for( y = this->m_YDim - 1 ; y >= 0 ; y-- )
	{
	  this->m_fstar_direction = 2 ;
	  for( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      this->VisitVoxel ( x, y, z ,3,i) ;
	    }
	  this->m_fstar_direction = 3 ;
	  for( x = this->m_XDim - 1 ; x >= 0 ; x-- )
	    {
	      this->VisitVoxel ( x, y, z ,4,i) ;
	    }
	}
    }

  std::cout << "Direction flip. " << this->m_NChanged << " voxels updated so far in this pass." << std::endl ;
  for( z = this->m_ZDim - 1 ; z >= 0 ; z-- )
    {
      for( y = 0 ; y < this->m_YDim ; y++ )
	{
	  this->m_fstar_direction = 4 ;
	  for( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      this->VisitVoxel ( x, y, z ,5,i) ;
	    }
	  this->m_fstar_direction = 5 ;
	  for( x = this->m_XDim - 1 ; x >= 0 ; x-- )
	    {
	      this->VisitVoxel ( x, y, z ,6,i) ;
	    }
	}
      for( y = this->m_YDim - 1 ; y >= 0 ; y-- )
	{
	  this->m_fstar_direction = 6 ;
	  for( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      this->VisitVoxel ( x, y, z ,7,i) ;
	    }
	  this->m_fstar_direction = 7 ;
	  for( x = this->m_XDim - 1 ; x >= 0 ; x-- )
	    {
	      this->VisitVoxel ( x, y, z ,8,i) ;
	    }
	}
    }
}

// locate all voxels that are in the source region, store their coordinates in a table
void ODFCost::FindSourceVoxels ()
{
  unsigned long index ;
  ODFCost::CoordinateType currentVoxel ;
  this->m_NumberOfSourceVoxels = 0 ;
  this->m_SourceVoxels.clear () ;

  for ( unsigned long z = 0 ; z < this->m_ZDim ; z++ )
    {
      for ( unsigned long y = 0 ; y < this->m_YDim ; y++ )
	{
	  for ( unsigned long x = 0 ; x < this->m_XDim ; x++ )
	    {
	       
	      index = x + y * this->m_XDim + z * this->m_Slice ;
	      if ( this->m_SourceArray[index] )
		{
		  currentVoxel[0] = x ;
		  currentVoxel[1] = y ;
		  currentVoxel[2] = z ;
		  this->m_SourceVoxels.push_back ( currentVoxel ) ;
		  this->m_NumberOfSourceVoxels ++ ;
		}
	    }
	}
    }
}

// allocate and initialize cost images
void ODFCost::CreateCostMaps () 
{
   //create empty output image
   this->m_CostImage = ODFImageType::New () ;
   this->m_AvgCostImage = ODFImageType::New () ;
   this->m_FinalCostImage = DoubleImageType::New () ;
   this->m_LengthImage = ODFImageType::New () ;
   this->m_OrigImage = ODFImageType::New () ;
   this->m_FinalOrigImage = DoubleImageType::New () ;
   this->m_FinalLengthImage = DoubleImageType::New () ;


   this->m_FinalPathImage = DoubleImageType::New () ;
   this->m_NeighborImage = ODFImageType::New () ;
   this->m_ODFVoxInImage = ODFImageType::New () ;
   this->m_ODFNeiOutImage = ODFImageType::New () ;
   this->m_PenInOutImage = ODFImageType::New () ;
   this->m_PenInDnImage = ODFImageType::New () ;
   this->m_PenOutDnImage = ODFImageType::New () ;

   this->m_FinalODFVoxInImage = DoubleImageType::New () ;
   this->m_FinalODFNeiOutImage = DoubleImageType::New () ;
   this->m_FinalPenInOutImage = DoubleImageType::New () ;
   this->m_FinalPenInDnImage = DoubleImageType::New () ;
   this->m_FinalPenOutDnImage = DoubleImageType::New () ;

   //added by Wenyu

   this->m_OrigInDirImage = DoubleImageType::New () ;
   this->m_OrigOutDirImage = DoubleImageType::New () ;

   itk::VariableLengthVector< ODFType > f ( this->m_NumberOfDirs ) ;
   for ( unsigned long i = 0 ; i < this->m_NumberOfDirs ; i++ ) 
     { 
       f[i] = INF; 
     }

   ODFImageType::IndexType start;
   ODFImageType::SizeType  size;
   start[0] = 0;   start[1] = 0;   start[2] = 0;
   size = this->m_ODFImage->GetLargestPossibleRegion().GetSize() ;

   this->m_XDim = size[0] ;
   this->m_YDim = size[1] ;
   this->m_ZDim = size[2] ;
   this->m_Slice = this->m_YDim * this->m_XDim ;

   ODFImageType::SpacingType spacing;
   spacing = this->m_ODFImage->GetSpacing () ;
   this->m_CostImage->SetSpacing ( spacing ) ;
   this->m_FinalCostImage->SetSpacing ( spacing ) ;
   this->m_LengthImage->SetSpacing ( spacing ) ;
   this->m_AvgCostImage->SetSpacing ( spacing ) ;
   this->m_OrigImage->SetSpacing ( spacing ) ;
   this->m_FinalOrigImage->SetSpacing ( spacing ) ;
   this->m_FinalLengthImage->SetSpacing ( spacing ) ;

   this->m_FinalPathImage->SetSpacing ( spacing ) ;
   this->m_NeighborImage->SetSpacing ( spacing ) ;

   this->m_ODFVoxInImage->SetSpacing ( spacing ) ;
   this->m_ODFNeiOutImage->SetSpacing ( spacing ) ;
   this->m_PenInOutImage->SetSpacing ( spacing ) ;
   this->m_PenInDnImage->SetSpacing ( spacing ) ;
   this->m_PenOutDnImage->SetSpacing ( spacing ) ;

   this->m_FinalODFVoxInImage->SetSpacing ( spacing ) ;
   this->m_FinalODFNeiOutImage->SetSpacing ( spacing ) ;
   this->m_FinalPenInOutImage->SetSpacing ( spacing ) ;
   this->m_FinalPenInDnImage->SetSpacing ( spacing ) ;
   this->m_FinalPenOutDnImage->SetSpacing ( spacing ) ;

   //added by Wenyu
   this->m_OrigInDirImage->SetSpacing ( spacing ) ;
   this->m_OrigOutDirImage->SetSpacing ( spacing ) ;

   ODFImageType::PointType origin ;
   origin = this->m_ODFImage->GetOrigin () ;
   this->m_CostImage->SetOrigin ( origin ) ;
   this->m_FinalCostImage->SetOrigin ( origin ) ;
   this->m_LengthImage->SetOrigin ( origin ) ;
   this->m_AvgCostImage->SetOrigin ( origin ) ;
   this->m_OrigImage->SetOrigin ( origin ) ;
   this->m_FinalOrigImage->SetOrigin ( origin ) ;
   this->m_FinalLengthImage->SetOrigin ( origin ) ;

   this->m_FinalPathImage->SetOrigin ( origin ) ;
   this->m_NeighborImage->SetOrigin ( origin ) ;

   this->m_ODFVoxInImage->SetOrigin ( origin ) ;
   this->m_ODFNeiOutImage->SetOrigin ( origin ) ;
   this->m_PenInOutImage->SetOrigin ( origin ) ;
   this->m_PenInDnImage->SetOrigin ( origin ) ;
   this->m_PenOutDnImage->SetOrigin ( origin ) ;

   this->m_FinalODFVoxInImage->SetOrigin ( origin ) ;
   this->m_FinalODFNeiOutImage->SetOrigin ( origin ) ;
   this->m_FinalPenInOutImage->SetOrigin ( origin ) ;
   this->m_FinalPenInDnImage->SetOrigin ( origin ) ;
   this->m_FinalPenOutDnImage->SetOrigin ( origin ) ;


   this->m_OrigInDirImage->SetOrigin ( origin ) ;
   this->m_OrigOutDirImage->SetOrigin ( origin ) ;

   ODFImageType::RegionType region;
   region.SetSize ( size ) ;
   region.SetIndex ( start ) ;
   this->m_CostImage->SetVectorLength ( this->m_NumberOfDirs ) ;
   this->m_CostImage->SetRegions ( region ) ;
   this->m_CostImage->Allocate () ;
   this->m_CostImage->FillBuffer ( f ) ;

   this->m_AvgCostImage->SetVectorLength ( this->m_NumberOfDirs ) ;
   this->m_AvgCostImage->SetRegions ( region ) ;
   this->m_AvgCostImage->Allocate () ;
   this->m_AvgCostImage->FillBuffer ( f ) ;

   this->m_LengthImage->SetVectorLength ( this->m_NumberOfDirs ) ;
   this->m_LengthImage->SetRegions ( region ) ;
   this->m_LengthImage->Allocate () ;
   this->m_LengthImage->FillBuffer ( f ) ;

   this->m_OrigImage->SetVectorLength ( this->m_NumberOfDirs ) ;
   this->m_OrigImage->SetRegions ( region ) ;
   this->m_OrigImage->Allocate () ;
   this->m_OrigImage->FillBuffer ( f ) ;

   this->m_FinalCostImage->SetRegions ( region ) ;
   this->m_FinalCostImage->Allocate () ;
   this->m_FinalCostImage->FillBuffer ( 0 ) ;

   this->m_FinalOrigImage->SetRegions ( region ) ;
   this->m_FinalOrigImage->Allocate () ;
   this->m_FinalOrigImage->FillBuffer ( 0 ) ;

   this->m_FinalLengthImage->SetRegions ( region ) ;
   this->m_FinalLengthImage->Allocate () ;
   this->m_FinalLengthImage->FillBuffer ( 0 ) ;

   this->m_FinalPathImage->SetRegions ( region ) ;
   this->m_FinalPathImage->Allocate () ;
   this->m_FinalPathImage->FillBuffer ( 0 ) ;

   this->m_NeighborImage->SetVectorLength ( this->m_NumberOfDirs ) ;
   this->m_NeighborImage->SetRegions ( region ) ;
   this->m_NeighborImage->Allocate () ;
   this->m_NeighborImage->FillBuffer ( f ) ;

   this->m_ODFVoxInImage->SetVectorLength ( this->m_NumberOfDirs ) ;
   this->m_ODFVoxInImage->SetRegions ( region ) ;
   this->m_ODFVoxInImage->Allocate () ;
   this->m_ODFVoxInImage->FillBuffer ( f ) ;

   this->m_ODFNeiOutImage->SetVectorLength ( this->m_NumberOfDirs ) ;
   this->m_ODFNeiOutImage->SetRegions ( region ) ;
   this->m_ODFNeiOutImage->Allocate () ;
   this->m_ODFNeiOutImage->FillBuffer ( f ) ;

   this->m_PenInOutImage->SetVectorLength ( this->m_NumberOfDirs ) ;
   this->m_PenInOutImage->SetRegions ( region ) ;
   this->m_PenInOutImage->Allocate () ;
   this->m_PenInOutImage->FillBuffer ( f ) ;

   this->m_PenInDnImage->SetVectorLength ( this->m_NumberOfDirs ) ;
   this->m_PenInDnImage->SetRegions ( region ) ;
   this->m_PenInDnImage->Allocate () ;
   this->m_PenInDnImage->FillBuffer ( f ) ;

   this->m_PenOutDnImage->SetVectorLength ( this->m_NumberOfDirs ) ;
   this->m_PenOutDnImage->SetRegions ( region ) ;
   this->m_PenOutDnImage->Allocate () ;
   this->m_PenOutDnImage->FillBuffer ( f ) ;

   this->m_FinalODFVoxInImage->SetRegions ( region ) ;
   this->m_FinalODFVoxInImage->Allocate () ;
   this->m_FinalODFVoxInImage->FillBuffer ( 0 ) ;

   this->m_FinalODFNeiOutImage->SetRegions ( region ) ;
   this->m_FinalODFNeiOutImage->Allocate () ;
   this->m_FinalODFNeiOutImage->FillBuffer ( 0 ) ;

   this->m_FinalPenInOutImage->SetRegions ( region ) ;
   this->m_FinalPenInOutImage->Allocate () ;
   this->m_FinalPenInOutImage->FillBuffer ( 0 ) ;

   this->m_FinalPenInDnImage->SetRegions ( region ) ;
   this->m_FinalPenInDnImage->Allocate () ;
   this->m_FinalPenInDnImage->FillBuffer ( 0 ) ;

   this->m_FinalPenOutDnImage->SetRegions ( region ) ;
   this->m_FinalPenOutDnImage->Allocate () ;
   this->m_FinalPenOutDnImage->FillBuffer ( 0 ) ;

   //added by Wenyu

   this->m_OrigInDirImage->SetRegions ( region ) ;
   this->m_OrigInDirImage->Allocate () ;
   this->m_OrigInDirImage->FillBuffer ( 0 ) ;


   this->m_OrigOutDirImage->SetRegions ( region ) ;
   this->m_OrigOutDirImage->Allocate () ;
   this->m_OrigOutDirImage->FillBuffer ( 0 ) ;

   this->m_CostArray = this->m_CostImage->GetPixelContainer()->GetBufferPointer();
   this->m_AvgCostArray = this->m_AvgCostImage->GetPixelContainer()->GetBufferPointer();
   this->m_FinalCostArray = this->m_FinalCostImage->GetPixelContainer()->GetBufferPointer() ;
   this->m_LengthArray = this->m_LengthImage->GetPixelContainer()->GetBufferPointer();
   this->m_OrigArray = this->m_OrigImage->GetPixelContainer()->GetBufferPointer();
   this->m_FinalOrigArray = this->m_FinalOrigImage->GetPixelContainer()->GetBufferPointer();
   this->m_FinalLengthArray = this->m_FinalLengthImage->GetPixelContainer()->GetBufferPointer();

   this->m_NeighborArray = this->m_NeighborImage->GetPixelContainer()->GetBufferPointer();

   this->m_FinalPathArray = this->m_FinalPathImage->GetPixelContainer()->GetBufferPointer();

   this->m_ODFVoxInArray = this->m_ODFVoxInImage->GetPixelContainer()->GetBufferPointer();
   this->m_ODFNeiOutArray = this->m_ODFNeiOutImage->GetPixelContainer()->GetBufferPointer();
   this->m_PenInOutArray = this->m_PenInOutImage->GetPixelContainer()->GetBufferPointer();
   this->m_PenInDnArray = this->m_PenInDnImage->GetPixelContainer()->GetBufferPointer();
   this->m_PenOutDnArray = this->m_PenOutDnImage->GetPixelContainer()->GetBufferPointer();

   this->m_FinalODFVoxInArray = this->m_FinalODFVoxInImage->GetPixelContainer()->GetBufferPointer();
   this->m_FinalODFNeiOutArray = this->m_FinalODFNeiOutImage->GetPixelContainer()->GetBufferPointer();
   this->m_FinalPenInOutArray = this->m_FinalPenInOutImage->GetPixelContainer()->GetBufferPointer();
   this->m_FinalPenInDnArray = this->m_FinalPenInDnImage->GetPixelContainer()->GetBufferPointer();
   this->m_FinalPenOutDnArray = this->m_FinalPenOutDnImage->GetPixelContainer()->GetBufferPointer();


 //added by Wenyu
   this->m_OrigInDir = this->m_OrigInDirImage->GetPixelContainer()->GetBufferPointer();
   this->m_OrigOutDir = this->m_OrigOutDirImage->GetPixelContainer()->GetBufferPointer();

   this->m_Cost = this->m_CostImage->GetPixelContainer()->GetBufferPointer();
}

// source voxels have 0 cost
void ODFCost::InitializeSourceVoxels () 
{
  ODFCost::CoordinateType currentVoxel ;
  unsigned long index, indexTimesNDirs ;
  unsigned int counter = 1 ;

  
  unsigned long int x, y, z;
  index = 0, indexTimesNDirs = 0 ;
  for ( z = 0 ; z < this->m_ZDim ; z++ )
    {
      for ( y = 0 ; y < this->m_YDim ; y++ )
	{
	  for ( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      index =  x + y * this->m_XDim + z * this->m_Slice ;
	      
	      this->m_OrigInDir[index] = this->m_NumberOfDirs;
	      this->m_OrigOutDir[index] = this->m_NumberOfDirs;
	      
	      indexTimesNDirs = index * this->m_NumberOfDirs ;
	      for ( unsigned long dir = 0 ; dir < this->m_NumberOfDirs ; dir++ )
		{
		  this->m_CostArray[indexTimesNDirs + dir] = INF;
		  this->m_AvgCostArray[indexTimesNDirs + dir] = INF;
		  this->m_LengthArray[indexTimesNDirs+dir] = 0   ;
		  this->m_OrigArray[indexTimesNDirs+dir] = this->m_NumberOfDirs ; 		
		}
	      
	      this->m_FinalCostArray[index] = INF ;
	    }
	}
    } 
  
  for ( unsigned long i = 0 ; i < this->m_NumberOfSourceVoxels ; i++ )
    {
      currentVoxel = this->m_SourceVoxels[i] ;
      
      index = currentVoxel[0] + (currentVoxel[1]) * this->m_XDim + (currentVoxel[2]) * this->m_Slice ;
      indexTimesNDirs = index * this->m_NumberOfDirs ;
      for ( unsigned long dir = 0 ; dir < this->m_NumberOfDirs ; dir ++ )
	{
	  this->m_CostArray[indexTimesNDirs+dir] = 0 ;
	  this->m_AvgCostArray[indexTimesNDirs+dir] = 0 ;
	  this->m_LengthArray[indexTimesNDirs+dir] = 0 ;
	  this->m_OrigArray[indexTimesNDirs+dir] = counter ;
	}
      this->m_FinalCostArray[index] = 0 ;
      counter++ ;
    }
}

// visit voxel, updating cost
void ODFCost::VisitVoxel ( unsigned long x, unsigned long y, unsigned long z ,int iter, int pass)
{
  unsigned long voxelIndex = x + y * this->m_XDim + z * this->m_Slice ;
  unsigned long voxelIndexTimesNDirs = voxelIndex * this->m_NumberOfDirs ;
//  if ( this->m_OrigFAArray[voxelIndex] < 0.05 ||  this->m_SourceArray[voxelIndex] )
  // If outside the White Matter
  if ( this->m_WMmaskArray[voxelIndex] == 0 || this->m_SourceArray[voxelIndex] ) // Added by Adrien Kaiser
    {
      // not worth the computation time TODO: add as command line argument => done by Adrien Kaiser 11/15/2012
      return ;
    }
//std::cout<<x<<","<<y<<","<<z<<std::endl;
  unsigned long neighborIndex, neighborIndexTimesNDirs ;
  unsigned short neighborIterator ;
  double newStepCost ;

  unsigned int origInDir ;
  unsigned int origOutDir ;

  // instead of all 26 neighbors, only visit the neighbors as indicated by fstar for our current traversal direction
  for ( unsigned short n = 0 ; n < this->m_fstar_updates[this->m_fstar_direction][0] ; n++ )
    {
      neighborIterator = this->m_fstar_updates[this->m_fstar_direction][n+1] ;

      if ( ! this->m_NeighborValidTable[voxelIndex][neighborIterator] ) continue ;
      neighborIndex = this->m_NeighborIndexTable[voxelIndex][neighborIterator] ;

      neighborIndexTimesNDirs = neighborIndex * this->m_NumberOfDirs ;
      double dn = this->m_DNTable[neighborIterator] ;

      origInDir = this->m_OrigInDir[ neighborIndex ];
      origOutDir = this->m_OrigOutDir[ neighborIndex];
   
	// for each incoming dir
	for ( unsigned int inDir = 0 ; inDir < this->m_NumberOfDirs ; inDir++ )  
	{
	  unsigned int voxelAndIn = voxelIndexTimesNDirs + inDir ;

	  // for each outgoing dir
	   for ( unsigned int outDir = 0 ; outDir < this->m_NumberOfDirs ; outDir++ )
	    {
	      unsigned int neighborAndOut = neighborIndexTimesNDirs+outDir ;

	      //I use this constraints for mouse data, it's unneccessary for synthetic data
	      if(this->m_DirectionTable[inDir][origInDir] && this->m_DirectionTable[outDir][origOutDir])
	     	{
 		  double preAngleCost =  this->m_AnglePenaltyTable[inDir][origInDir] + this->m_AnglePenaltyTable[outDir][origOutDir] ;

		  double f_odf1 =  this->m_FinslerTable[voxelAndIn];  
		  double f_odf2 = this->m_FinslerTable[neighborAndOut];
		  double f_odf = f_odf1 + f_odf2;

		  //For synthetic data, we do not need the preAngleCost part
		 double cost = this->m_NeighborAnglePenaltyTable[inDir][neighborIterator] + this->m_NeighborAnglePenaltyTable[outDir][neighborIterator]+this->m_AnglePenaltyTable[inDir][outDir]+preAngleCost;
		  //double cost = this->m_NeighborAnglePenaltyTable[inDir][neighborIterator] + this->m_NeighborAnglePenaltyTable[outDir][neighborIterator]+this->m_AnglePenaltyTable[inDir][outDir];

		  newStepCost = (  this->m_alpha*f_odf + cost ) * dn ;

	// l.651: void ODFCost::UpdateCost ( unsigned long index, unsigned long voxelIndex, double newCost,           double oldAvgCost,                   double newLength,                    unsigned int newOrig,       unsigned int indir, unsigned int outdir, unsigned int neighbor, unsigned int neighborAndOut)
		  if(this->UpdateCost ( voxelAndIn, voxelIndex, newStepCost + this->m_CostArray[neighborAndOut], this->m_AvgCostArray[neighborAndOut], dn+this->m_LengthArray[neighborAndOut], this->m_OrigArray[neighborAndOut], inDir, outDir, neighborIterator ,neighborAndOut))
			{
			   this->m_ODFVoxInArray[voxelAndIn] = f_odf1;
			   this->m_ODFNeiOutArray[voxelAndIn] = f_odf2;
			   this->m_PenInOutArray[voxelAndIn] = this->m_AnglePenaltyTable[inDir][outDir];
			   this->m_PenInDnArray[voxelAndIn] = this->m_NeighborAnglePenaltyTable[inDir][neighborIterator];
			   this->m_PenOutDnArray[voxelAndIn] = this->m_NeighborAnglePenaltyTable[outDir][neighborIterator];
			}

		} // if mouse data
	    }// for each outgoing dir
	} // for each incoming dir
    } // for each neighbor

//  std::cout<<"["<<x<<","<<y<<","<<z<<"] (iter="<<iter<<", inDir="<<inDir <<") AvgCost( "<<voxelAndIn<<" )= "<< this->m_AvgCostArray[voxelAndIn] << std::endl;

if( x==46 && y==47 && z==38 ) // 2 source= [46,48,38] => get values in [46,47,38] -> minDir = 12 and min =0.17414
{
    unsigned int voxelAndIn = voxelIndexTimesNDirs + 12 ;
    std::cout<<pass<<","<<iter<<","<< this->m_AvgCostArray[voxelAndIn] << ","<< this->m_CostArray[voxelAndIn] <<","<< this->m_LengthArray[voxelAndIn] <<","<< this->m_OrigArray[voxelAndIn] <<std::endl;
}

if( x==46 && y==48 && z==38 ) // 1 source= [46,47,38] => get values in [46,48,38] -> minDir = 1 and min =0.17004
{
	double minCost=INF;
	double minDir;
	for ( unsigned int inDir = 0 ; inDir < this->m_NumberOfDirs ; inDir++ )  
	{
	  unsigned int voxelAndIn = voxelIndexTimesNDirs + inDir ;
	   if( this->m_AvgCostArray[voxelAndIn] < minCost) minDir=inDir;
	}

    unsigned int voxelAndIn = voxelIndexTimesNDirs + minDir ;
    std::cout<<pass<<","<<iter<<","<<minDir<<","<< this->m_AvgCostArray[voxelAndIn] << ","<< this->m_CostArray[voxelAndIn] <<","<< this->m_LengthArray[voxelAndIn] <<","<< this->m_OrigArray[voxelAndIn] <<std::endl;
}

} // end function

void ODFCost::PrecomputeDirection()
{
  std::cout<<"Compute direction table"<<std::endl;
  this->m_DirectionTable.resize ( this->m_NumberOfDirs + 1) ;
  for(unsigned int i = 0; i < this->m_NumberOfDirs ; i ++)
    {
      this->m_DirectionTable[i].resize ( this->m_NumberOfDirs + 1) ;
      for(unsigned int j=0; j < this->m_NumberOfDirs; j++)
	{
	  ODFCost::CoordinateType jdir = this->m_CoordinateTable[j];
	  ODFCost::CoordinateType idir = this->m_CoordinateTable[i];  
	  double magi = sqrt(idir[0]*idir[0]+idir[1]*idir[1]+idir[2]*idir[2]);
	  double magj = sqrt(jdir[0]*jdir[0]+jdir[1]*jdir[1]+jdir[2]*jdir[2]);
	  double sum = idir[0]*jdir[0]+idir[1]*jdir[1]+idir[2]*jdir[2];
	  double angle = acos(sum/(magi*magj));
	  if( angle < 0.5*PI )
	    this->m_DirectionTable[i][j] = 1;
	  else
	    this->m_DirectionTable[i][j] = 0; 
	}
      this->m_DirectionTable[i][this->m_NumberOfDirs] = 1;
    }
  this->m_DirectionTable[this->m_NumberOfDirs].resize ( this->m_NumberOfDirs + 1) ;
  for(unsigned int k = 0; k <= this->m_NumberOfDirs; k ++ )
    {
      this->m_DirectionTable[this->m_NumberOfDirs][k] = 1;
    }
}

//added by Wenyu
void ODFCost::CreateFinslerTable()
{
  std::cout<<"Compute Finsler Table"<<std::endl;
  unsigned long int x, y, z, index = 0, indexTimesNDirs = 0 ;
  ODFCost::CoordinateType idir, jdir ;
  unsigned long curdir; 
  double f_ODF;

  double min_fodf = 1, max_fodf = 0 ;
  index = -1 ;

  indexTimesNDirs = -this->m_NumberOfDirs ;
  for ( z = 0 ; z < this->m_ZDim ; z++ )
    {
      for ( y = 0 ; y < this->m_YDim ; y++ )
	{
	  for ( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      index++ ;
	      
	      indexTimesNDirs += this->m_NumberOfDirs ;
	      for ( curdir = 0 ; curdir < this->m_NumberOfDirs ; curdir++ )
		{
		  idir = this->m_CoordinateTable[curdir] ;		  	  

		  double curODF = this->m_ODFArray [ indexTimesNDirs+curdir ] ; 
		  if ( curODF != 0 ) 
		    {
		      double sumODF_proj = 0 ;
		      curODF *= 0.25*this->m_NumberOfDirs; // NOT SURE WHAT THE APPROPRIATE COEFFICIENT HERE SHOULD REALLY BE
		      for ( unsigned long dir = 0 ; dir < this->m_NumberOfDirs ; dir++ )
			{
			  jdir = this->m_CoordinateTable [dir] ;
			  double dirODF = this->m_ODFArray [indexTimesNDirs+dir];
			  double cosAlpha = idir[0]*jdir[0]+idir[1]*jdir[1]+idir[2]*jdir[2] ;
			  double sinAlpha2 = 1 - cosAlpha * cosAlpha ;
			  double sinAlpha = sinAlpha2 < 0 ? 0 : sqrt ( sinAlpha2 ) ;
			  sumODF_proj += dirODF * sinAlpha ;
			}

		      double finsler = curODF/sumODF_proj ;
		      if ( finsler < 0 ) 
			{
			  std::cout << "Something is very wrong. This should never be negative. Check your ODF file. " << std::endl ;
			  exit ( 0 ) ;
			}
		      else if ( finsler > 1 ) 
			{
			  std::cout << curODF << " " << sumODF_proj << std::endl ;
			  std::cout << "Something is very wrong. This should never be higher than 1. Check your ODF file and/or the normalization. " << std::endl ;
			  exit ( 0 ) ;
			}
		      f_ODF = 1 - finsler ;

		    }
		  else // curodf== 0
		    {
		      f_ODF = 1 ;
		    }

		  if ( min_fodf > f_ODF )
		    min_fodf = f_ODF ; 
		  if ( max_fodf < f_ODF )
		    max_fodf = f_ODF ;
		  f_ODF = f_ODF*f_ODF*f_ODF*f_ODF*f_ODF*f_ODF;
		  this->m_FinslerTable.push_back ( f_ODF ) ;
		}
	    }
	}
    }
  
  std::cout<<"min_fodf is "<<min_fodf<<","<<"max_fodf is "<<max_fodf<<std::endl;
}

// for each voxel, the final cost is the minimum cost among the costs associated with all the directions
void ODFCost::ComputeFinalCosts ()
{

// 106 * 106 * 76 = 853936
unsigned long FinalNeighbors[853936];

  unsigned long int x, y, z, index = 0, indexTimesNDirs = 0 ;
  double currentCost, minCost, minLength ;
  unsigned int minOrig;
std::cout<< "this->m_XDim= "<<this->m_XDim<<" | this->m_Slice= "<<this->m_Slice<<" | this->m_NumberOfDirs= "<<this->m_NumberOfDirs<<std::endl;

  for ( z = 0 ; z < this->m_ZDim ; z++ )
    {
      for ( y = 0 ; y < this->m_YDim ; y++ )
	{
	  for ( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      minCost = INF ;
	      minOrig = -1 ;
	      minLength = -1 ;
	      unsigned long mindir = 0;
	      double minODF = 2;
	      double odfcost;
	      unsigned int odfdir;

              double minODFVoxIn=0;
              double minODFNeiOut=0;
              double minPenInOut=0;
              double minPenInDn=0;
              double minPenOutDn=0;

	      double sum = 0;
	      for ( unsigned long dir = 0 ; dir < this->m_NumberOfDirs ; dir++ )
		{
		 
/*
		if(x==46 && y==48 && z==38) // to display on sphere
		{
		  std::cout<<dir<<","<<this->m_ODFVoxInArray[indexTimesNDirs+dir]<<","<<this->m_ODFNeiOutArray[indexTimesNDirs+dir]<<std::endl;
		}
*/
		  currentCost = this->m_AvgCostArray[indexTimesNDirs+dir];
		  sum += currentCost;		    

		  if ( ( currentCost < minCost ) && ( currentCost >= 0 ) )
		    {
		      mindir = dir;
		      minCost = currentCost ;
		      minOrig = this->m_OrigArray[indexTimesNDirs+dir] ;
		      minLength = this->m_LengthArray[indexTimesNDirs+dir] ;

		  minODFVoxIn=this->m_ODFVoxInArray[indexTimesNDirs+dir] ;
		  minODFNeiOut=this->m_ODFNeiOutArray[indexTimesNDirs+dir] ;
		  minPenInOut=this->m_PenInOutArray[indexTimesNDirs+dir] ;
		  minPenInDn=this->m_PenInDnArray[indexTimesNDirs+dir] ;
		  minPenOutDn=this->m_PenOutDnArray[indexTimesNDirs+dir] ;
	    }

		  if(minODF > this->m_FinslerTable[indexTimesNDirs+ dir])
		    {
		      minODF =  this->m_FinslerTable[indexTimesNDirs+ dir];
		      odfdir = dir;
		      odfcost = currentCost;
		    }
		  
		} 

	      if ( minCost == INF ) minCost = -1 ;
	     
	      this->m_FinalCostArray[index] = minCost ;
	      
	      this->m_FinalOrigArray[index] = minOrig ;
	      this->m_FinalLengthArray[index] = minLength ;

		if(x==46 && y==48 && z==38) // to display on sphere
		{
		  std::cout<<"46,48,38: final length="<<minLength<<std::endl;
		}

		if(minODFVoxIn!=INF)   this->m_FinalODFVoxInArray[index] = minODFVoxIn ;
		if(minODFNeiOut!=INF)  this->m_FinalODFNeiOutArray[index] = minODFNeiOut ;
		if(minPenInOut!=INF)  this->m_FinalPenInOutArray[index] = minPenInOut ;
		if(minPenInDn!=INF)  this->m_FinalPenInDnArray[index] = minPenInDn ;
		if(minPenOutDn!=INF)  this->m_FinalPenOutDnArray[index] = minPenOutDn ;

	      index++ ;
	      indexTimesNDirs += this->m_NumberOfDirs ;

		  unsigned long voxelIndex = x + y * this->m_XDim + z * this->m_Slice ;
		  unsigned long voxelIndexTimesNDirs = voxelIndex * this->m_NumberOfDirs ;
		  unsigned int  voxelAndIn = voxelIndexTimesNDirs + mindir ;

		unsigned long neighborWithDirs = this->m_NeighborArray[voxelAndIn]; // neighbor with directions
		unsigned long neighborWithNODirs = (int)(neighborWithDirs/this->m_NumberOfDirs) ; // neighbor with NO directions : ex : voxelAndIn -> voxelIndex // r=a-b*q // a=154, b=5 : q=30, r=4 : neighborWithNODirs = q (= voxelIndex)

		FinalNeighbors[voxelIndex]=neighborWithNODirs;


		if( (x==46 && y==47 && z==38) || (x==46 && y==48 && z==38) ) // (pair 1/2)
		{
		std::cout<<"Final Cost= "<< minCost << std::endl;
		}

	    }
	}
    }
	unsigned long int xs=46;
	unsigned long int ys=47;
	unsigned long int zs=38;
	unsigned long int xt=46;
	unsigned long int yt=48;
	unsigned long int zt=38;

	unsigned long source = xs + ys * this->m_XDim + zs * this->m_Slice ;
	unsigned long target = xt + yt * this->m_XDim + zt * this->m_Slice ;

	unsigned long neighbor= target;
	int zn=(int)(neighbor/this->m_Slice);
	int yn=(int)( (neighbor-zn*this->m_Slice) / this->m_XDim );
	int xn= neighbor - zn*this->m_Slice - yn*this->m_XDim;
	std::cout<<neighbor<<" \t| "<<xn<<","<<yn<<","<<zn<<std::endl;
	int Nbneighbors=0;
	while (neighbor!=source && Nbneighbors<1000000)
	{
	m_FinalPathArray[neighbor]=1;
	neighbor=FinalNeighbors[neighbor];
	std::cout<<" v "<<std::endl;
	zn=(int)(neighbor/this->m_Slice);
	yn=(int)( (neighbor-zn*this->m_Slice) / this->m_XDim );
	xn= neighbor - zn*this->m_Slice - yn*this->m_XDim;
	std::cout<<neighbor<<" \t| "<<xn<<","<<yn<<","<<zn<<std::endl;
	Nbneighbors++;
	}
	std::cout<<"Nbneighbors="<<Nbneighbors<<std::endl;

}

// get the neighbor index (x+y*xDim+z*xDim*yDim) given current voxel and 0 <= neighbor < 27
// if the neighbor is out of the image volume, "valid" becomes false
unsigned long ODFCost::GetNeighbor ( unsigned long x, unsigned long y, unsigned long z, unsigned short neighbor, bool &valid) 
{
  short nx, ny, nz ;
  nx = neighbor % 3;
  neighbor = ( neighbor - nx) / 3;
  nx-- ;
  ny = neighbor % 3 ;
  nz = ( neighbor - ny ) / 3 - 1 ;
  ny-- ;
  
  short neighborX = nx + x ;
  short neighborY = ny + y ;
  short neighborZ = nz + z ;
  
  valid = true ;
  if ( neighborX < 0 || neighborX >= this->m_XDim || neighborY < 0 || neighborY >= this->m_YDim || neighborZ < 0 || neighborZ >= this->m_ZDim )
    valid = false ;

  return neighborX + neighborY * this->m_XDim + neighborZ * this->m_Slice ;
}


void ODFCost::PrecomputeNeighbors ()
{
  std::cout << "Precompute neighbors" << std::endl ;
  unsigned long nVoxels = this->m_XDim * this->m_YDim * this->m_ZDim ;   
  this->m_NeighborIndexTable.resize ( nVoxels ) ;
  this->m_NeighborValidTable.resize ( nVoxels ) ;

  unsigned long neighborIndex ;
  bool neighborValid ;

  for ( unsigned int z = 0 ; z < this->m_ZDim ; z++ )
    {
      for ( unsigned int y = 0 ; y < this->m_YDim ; y++ )
	{
	  for ( unsigned int x = 0 ; x < this->m_XDim ; x++ )
	    {
	      unsigned long currentVoxel = x + y * this->m_XDim + z * this->m_Slice ;
	      this->m_NeighborIndexTable[currentVoxel].resize ( 27 ) ;
	      this->m_NeighborValidTable[currentVoxel].resize ( 27 ) ;

	      for ( unsigned int neighbor = 0 ; neighbor < 27 ; neighbor++ )
		{
		  neighborIndex = this->GetNeighbor ( x, y, z, neighbor, neighborValid ) ;
		  this->m_NeighborIndexTable[currentVoxel][neighbor] = neighborIndex ;
		  this->m_NeighborValidTable[currentVoxel][neighbor] = neighborValid ;
		}
	    }
	}
   }
}

// if the current path is cheaper than the current best path, update
//modified by Wenyu
//void ODFCost::UpdateCost ( unsigned long index, double newCost, double newLength, unsigned int newOrig, unsigned int dir)
bool ODFCost::UpdateCost ( unsigned long index, unsigned long voxelIndex, double newCost,double oldAvgCost, double newLength, unsigned int newOrig, unsigned int indir,unsigned int outdir, unsigned int neighbor, unsigned int neighborAndOut)
{ // index= VoxelAndIn
  double currentAvgCost = this->m_AvgCostArray[index] ;
  double newAvgCost = newCost / newLength ;
  
  if ( newAvgCost < oldAvgCost)
    {
      newAvgCost = oldAvgCost;
    } 
  
  if ( newAvgCost*1.05 < currentAvgCost ) 
    {
     this->m_CostArray[index] = newCost ;
      this->m_Converged = false ;
      this->m_NChanged++ ;
      this->m_LengthArray[index] = newLength ;
      this->m_AvgCostArray[index] = newAvgCost ;
      this->m_OrigArray[index] = newOrig;
//std::cout<<"PONEY6"<<std::endl;
      // added by Adrien Kaiser
//      std::cout<<"voxel="<< index<<" \t| new neighbor="<<neighborAndOut<<std::endl;
      this->m_NeighborArray[index] = neighborAndOut; // index is VoxelAndIn
//std::cout<<"PONEY8"<<std::endl;
      //added by Wenyu
      this->m_OrigInDir[voxelIndex] = indir;
      this->m_OrigOutDir[voxelIndex] = outdir;

       return true;
    }
return false;
}

// the penalty associated with the angle discrepancy between two directions (arguments are indices from 0 to nDirs)
double ODFCost::AnglePenalty ( unsigned u, unsigned v )
{
  ODFCost::CoordinateType idir, jdir ;
  idir = this->m_CoordinateTable[u] ;
  jdir = this->m_CoordinateTable[v] ;
  
  double sum = idir[0] * jdir[0] + idir[1] * jdir[1] + idir[2] * jdir[2];
//if(sum<0) sum = -1*sum;
/*  double angle = acos(sum); // 0 < angle < PI/2
  if ( angle>PI/2 ) return 1.;
  else if ( angle<PI/2 ) return 0.;
  else return 0.5; // if ( angle==PI/2 )
*/

  double e = 0.5 * ( 1. - sum ) ;
 
  return e;
}

// penalty associated with the angle discrepancy between two directions ( 0 <= u < nDirs, 0 <= neighbor < 27 )
double ODFCost::NeighborAnglePenalty ( unsigned u, unsigned neighbor )
{
  ODFCost::CoordinateType idir, jdir, kdir ;
  idir = this->m_CoordinateTable[u] ;

  //std::cout<<"coornate of dir "<<u<<" : "<<idir[0]<<" "<<idir[1]<<" "<<idir[2]<<std::endl;

  short nx, ny, nz ;
  nx = neighbor % 3  ;
  short neighborTemp = ( neighbor - nx ) / 3;
  nx-- ;
  ny = neighborTemp % 3 ;
  nz = ( neighborTemp - ny ) / 3 - 1 ;
  ny-- ;
  
  jdir[0] = nx / this->m_DNTable[neighbor] ;
  jdir[1] = ny / this->m_DNTable[neighbor] ;
  jdir[2] = nz / this->m_DNTable[neighbor] ;
  
  double sum = idir[0] * jdir[0] + idir[1] * jdir[1] + idir[2] * jdir[2] ;
//  if(sum<0) sum = -1*sum;
/*
  double angle = acos(sum); // 0 < angle < PI/2
  if ( angle>PI/2 ) return 1.;
  else if ( angle<PI/2 ) return 0.;
  else return 0.5; // if ( angle==PI/2 )
*/
  double e = 0.5 * ( 1. - sum ) ;
  return e ;
  
}

double ODFCost::NeighborPenalty ( unsigned neighbor1, unsigned neighbor2 )
{
  ODFCost::CoordinateType idir, jdir ;

  short nx1, ny1, nz1 ;
  nx1 = neighbor1 % 3  ;
  short neighborTemp1 = ( neighbor1 - nx1 ) / 3;
  nx1-- ;
  ny1 = neighborTemp1 % 3 ;
  nz1 = ( neighborTemp1 - ny1 ) / 3 - 1 ;
  ny1-- ;
  
  idir[0] = nx1 / this->m_DNTable[neighbor1] ;
  idir[1] = ny1 / this->m_DNTable[neighbor1] ;
  idir[2] = nz1 / this->m_DNTable[neighbor1] ;

  short nx2, ny2, nz2 ;
  nx2 = neighbor2 % 3 ;
  short neighborTemp2 = ( neighbor2 - nx2 ) / 3;
  nx2-- ;
  ny2 = neighborTemp2 % 3 ;
  nz2 = ( neighborTemp2 - ny2 ) / 3 - 1 ;
  ny2-- ;
  
  jdir[0] = nx2 / this->m_DNTable[neighbor2] ;
  jdir[1] = ny2 / this->m_DNTable[neighbor2] ;
  jdir[2] = nz2 / this->m_DNTable[neighbor2] ;
 
  double sum = idir[0] * jdir[0] + idir[1] * jdir[1] + idir[2] * jdir[2] ;
  /*
  double angle = acos(sum); // 0 < angle < PI/2
  if ( angle>PI/2 ) return 1.;
  else if ( angle<PI/2 ) return 0.;
  else return 0.5; // if ( angle==PI/2 )
*/
  double e = 0.5 * ( 1. - sum ) ;
  return e ;
 
}

// a face-neighbor has a step length of 1, edge-neighbor \sqrt(2), vertex-neighbor \sqrt(3)
// 0 <= neighbor < 27
double ODFCost::stepLength ( unsigned neighbor ) 
{
  short nx, ny, nz ;
  nx = neighbor % 3  ;
  neighbor = ( neighbor - nx ) / 3;
  nx-- ;
  ny = neighbor % 3 ;
  nz = ( neighbor - ny ) / 3 - 1 ;
  ny-- ;

  short dn = abs ( nx ) + abs ( ny ) + abs ( nz ) ;
  if ( dn == 0 ) return 0 ;
  if ( dn == 1 ) return 1 ;
  if ( dn == 2 ) return sqrt ( 2 ) ;
  return sqrt ( 3 ) ;
}

// compute step lengths and store for faster access
void ODFCost::PrecomputeStepLength () 
{
  this->m_DNTable.resize ( 27 ) ;
  for ( unsigned i = 0 ; i < 27 ; i++ )
    {    
      this->m_DNTable[i] = this->stepLength ( i ) ;
    }
}


// compute angle penalties and store for faster access
void ODFCost::PrecomputeAnglePenalty () 
{
  std::cout<<"precompute angle penalty"<<std::endl;
  this->m_AnglePenaltyTable.resize ( this->m_NumberOfDirs+1 ) ;
  for ( unsigned i = 0 ; i <= this->m_NumberOfDirs ; i++ )
    {
      this->m_AnglePenaltyTable[i].resize ( this->m_NumberOfDirs + 1 ) ;
      for ( unsigned j = 0 ; j <= this->m_NumberOfDirs ; j++ )
	{
	  if(i ==  this->m_NumberOfDirs || j ==  this->m_NumberOfDirs)
	    {
	      this->m_AnglePenaltyTable[i][j] = 0;
	    }
	  else
	    {
	      this->m_AnglePenaltyTable[i][j] = this->AnglePenalty ( i, j ) ;
	    }
	 
	}
    }

  this->m_NeighborAnglePenaltyTable.resize ( this->m_NumberOfDirs+1 ) ;
  
  for ( unsigned i = 0 ; i <= this->m_NumberOfDirs ; i++ )
    {
      this->m_NeighborAnglePenaltyTable[i].resize ( 27 ) ;
      for ( unsigned j = 0 ; j < 27 ; j++ )
	{
	  if(i ==  this->m_NumberOfDirs)
	    {
	      this->m_NeighborAnglePenaltyTable[i][j] = 0 ;
	    }
	  else
	    {
	      this->m_NeighborAnglePenaltyTable[i][j] = this->NeighborAnglePenalty ( i, j ) ;
	      if(j == 13)
		{
		  this->m_NeighborAnglePenaltyTable[i][j] = 0;
		}
	     }
	}
      
    }
  this->m_NeighborPenaltyTable.resize ( 27 ) ;
  
  for ( unsigned i = 0 ; i < 27 ; i++ )
    {
      
      this->m_NeighborPenaltyTable[i].resize ( 27 ) ;
      for ( unsigned j = 0 ; j < 27 ; j++ )
	{
		 
	  this->m_NeighborPenaltyTable[i][j] = this->NeighborPenalty ( i, j ) ;
	  if(i==13 || j == 13)
	    {
	      this->m_NeighborPenaltyTable[i][j] = 0;
	    }
	  
	}
      
    }
}

// we want the odf values to be 0..1
void ODFCost::NormalizeODF () 
{
   double min = INF, max = -INF, current ;
   for ( unsigned long index = 0 ; index < this->m_ZDim * this->m_YDim * this->m_XDim * this->m_NumberOfDirs ; index++ )
     {
       current = this->m_ODFArray[index] ;
       if ( current < min ) min = current ;
       if ( current > max ) max = current ;
     }
   std::cout << "ODF min-max: " << min << " " << max << std::endl ;

  double sum ;
  unsigned long index, indexTimesNDir = 0 ;
  for ( index = 0 ; index < this->m_ZDim * this->m_YDim * this->m_XDim ; index++ )
    {
      sum = 0 ;
      for ( unsigned dir = 0 ; dir < this->m_NumberOfDirs ; dir++ )
	{
	  sum+= this->m_ODFArray[indexTimesNDir+dir] - min ;
	}
      for ( unsigned dir = 0 ; dir < this->m_NumberOfDirs ; dir++ )
	{
	  current = this->m_ODFArray[indexTimesNDir+dir] - min ;
	  current /= sum ;
	  this->m_ODFArray[indexTimesNDir+dir] = current ;
	}
      indexTimesNDir += this->m_NumberOfDirs ;
    }
}

// initialize the neighbors that need visiting along each fstar traversal direction
// this table is structured as follows:
// m_fstar_updates[i][0] = number (n) of neighbors that need visiting if we are on the i-direction (0<=i<8)
// m_fstar_updates[i][1..n] = the neighbors that need visiting
// the neighbors are indexed 0..27 (so 13 would be the voxel itself)
void ODFCost::FstarTraversalUpdates ()
{
  this->m_fstar_updates[0][0] = 13 ;
  this->m_fstar_updates[0][1] = 0 ;
  this->m_fstar_updates[0][2] = 1 ;
  this->m_fstar_updates[0][3] = 2 ;
  this->m_fstar_updates[0][4] = 3 ;
  this->m_fstar_updates[0][5] = 4 ;
  this->m_fstar_updates[0][6] = 5 ;
  this->m_fstar_updates[0][7] = 6 ;
  this->m_fstar_updates[0][8] = 7 ;
  this->m_fstar_updates[0][9] = 8 ;
  this->m_fstar_updates[0][10] = 9 ;
  this->m_fstar_updates[0][11] = 10 ;
  this->m_fstar_updates[0][12] = 11 ;
  this->m_fstar_updates[0][13] = 12 ;

  this->m_fstar_updates[1][0] = 1 ;
  this->m_fstar_updates[1][1] = 14 ;

  this->m_fstar_updates[2][0] = 4 ;
  this->m_fstar_updates[2][1] = 12 ;
  this->m_fstar_updates[2][2] = 15 ;
  this->m_fstar_updates[2][3] = 16 ;
  this->m_fstar_updates[2][4] = 17 ;

  this->m_fstar_updates[3][0] = 1 ;
  this->m_fstar_updates[3][1] = 14 ;

  this->m_fstar_updates[4][0] = 13 ;
  this->m_fstar_updates[4][1] = 14 ;
  this->m_fstar_updates[4][2] = 15 ;
  this->m_fstar_updates[4][3] = 16 ;
  this->m_fstar_updates[4][4] = 17 ;
  this->m_fstar_updates[4][5] = 18 ;
  this->m_fstar_updates[4][6] = 19 ;
  this->m_fstar_updates[4][7] = 20 ;
  this->m_fstar_updates[4][8] = 21 ;
  this->m_fstar_updates[4][9] = 22 ;
  this->m_fstar_updates[4][10] = 23 ;
  this->m_fstar_updates[4][11] = 24 ;
  this->m_fstar_updates[4][12] = 25 ;
  this->m_fstar_updates[4][13] = 26 ;

  this->m_fstar_updates[5][0] = 1 ;
  this->m_fstar_updates[5][1] = 14 ;

  this->m_fstar_updates[6][0] = 4 ;
  this->m_fstar_updates[6][1] = 12 ;
  this->m_fstar_updates[6][2] = 15 ;
  this->m_fstar_updates[6][3] = 16 ;
  this->m_fstar_updates[6][4] = 17 ;

  this->m_fstar_updates[7][0] = 1 ;
  this->m_fstar_updates[7][1] = 14 ;
}

