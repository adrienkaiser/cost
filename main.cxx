#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>

#include <iostream>
#include <string>
#include <fstream> // to open a file
#include "CostCLP.h"
#include "ODFReconstructor.h"
#include "ODFStreamline.h"
#include "ODFCost.h"
#include "define.h" // contains definitions for the tensor GK data file -> configured when ccmake

//ODF image definition
typedef double 					CostType;
typedef itk::VectorImage < CostType,3 > 	ODFImageType;
typedef itk::ImageFileReader < ODFImageType > 	ODFReaderType;
typedef itk::ImageFileWriter < ODFImageType > 	ODFWriterType;

typedef ODFReconstructor                        ODFReconstructorType ;
typedef ODFStreamline                           ODFStreamlineType ;
typedef ODFCost                                 ODFCostType ;

//Label image definition
typedef unsigned long                           ImageElementType;
typedef itk::Image < ImageElementType , 3 > 	ImageType;
typedef itk::ImageFileReader < ImageType > 	ReaderType;
typedef itk::ImageFileWriter < ImageType > 	WriterType;

//Double image definition
typedef itk::Image < double, 3 >                DoubleImageType ;
typedef itk::ImageFileReader < DoubleImageType > DoubleReaderType ;
typedef itk::ImageFileWriter < DoubleImageType > DoubleWriterType ;

//Color image definition
typedef itk::RGBPixel < unsigned char >        RGBPixelType ;
typedef itk::Image < RGBPixelType, 3 >         ColorImageType ;
typedef itk::ImageFileWriter < ColorImageType > ColorWriterType ;

std::string CreateDataFile() // Added by Adrien Kaiser : returns path to the new data file
{
  std::string Text = GK_Table_Text; // Defined in "define.h"

  // Write data file
  const char * HomeDir = getenv("HOME");
  std::string DataFile = (std::string)HomeDir + "/tensorGK.dat";
  if( access(HomeDir, W_OK) != 0 ) DataFile="tensorGK.dat";

  std::ofstream DataFileStream ( DataFile.c_str(), std::ios::out | std::ios::trunc);
  if (! DataFileStream) // error while opening
    {
      std::cout<<"Error creating the Data file \'"<< DataFile <<"\'. Abort."<<std::endl;
      exit(EXIT_FAILURE);
    }

  DataFileStream << Text << std::endl;

  DataFileStream.close();

  return DataFile;
}

int main(int argc, char * argv[])
{
  PARSE_ARGS ;

  //ODF image pointers
  ODFImageType::Pointer ODFCoefsImage 	= ODFImageType::New();
  ODFImageType::Pointer ODFImage 	= ODFImageType::New();
  ODFReaderType::Pointer ODFReader      = ODFReaderType::New();
  
  //added by Wenyu
  DoubleImageType::Pointer OrigFAImage            = DoubleImageType::New();
  DoubleReaderType::Pointer OrigFAReader          = DoubleReaderType::New();

  // WM mask  (Added by Adrien Kaiser)
  ImageType::Pointer WMmaskImage            = ImageType::New();
  ReaderType::Pointer WMmaskReader         = ReaderType::New();

  //Source (label) image pointers
  ImageType::Pointer SourceImage 	= ImageType::New();
  ReaderType::Pointer SourceReader      = ReaderType::New();

  // read source (seed) image	
  SourceReader->SetFileName ( source_image ) ;

  try
    {
      SourceReader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem reading the input source file: " << source_image << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
  SourceImage = SourceReader->GetOutput();

  //read original FA image
  OrigFAReader->SetFileName( fa_image );
  try
    {
      OrigFAReader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem reading the input fa file: " << fa_image << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
  OrigFAImage = OrigFAReader->GetOutput();

  //read WM mask (Added by Adrien Kaiser)
  WMmaskReader->SetFileName( wmmask_image );
  try
    {
      WMmaskReader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem reading the input wm mask file: " << wmmask_image << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
  WMmaskImage = WMmaskReader->GetOutput();

  // read odf image
  ODFReader->SetFileName ( odf_image ) ;
  try
    {
      ODFReader->Update() ; 
    }
  catch( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem reading the input ODF file: " << odf_image <<  std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
  ODFCoefsImage = ODFReader->GetOutput();

  //make the ODF usable
// Added by Adrien Kaiser : find the data GK file
  std::string datFile (GK_Table) ; // Added by Adrien Kaiser
  if ( access(datFile.c_str(), R_OK) != 0 ) // If the file is not found in the original location
    {
       std::cout<<"Gordon Kindlman's tensor table (file 'tensorGK.dat') not found: Please copy the file in the COST executable directory"<<std::endl;
/*      datFile=(std::string)getenv("HOME") + "/tensorGK.dat"; // found other location : home
      if ( access(datFile.c_str(), R_OK) != 0 ) // if the file is not found in other locations (home)
        {
          std::cout<<"Gordon Kindlman's tensor table (file 'tensorGK.dat') not found: A file will be created in your home directory for the following number of sampled on hemisphere: 50,100."<<std::endl;
          if ( numberOfSamplesOnHemisphere!=50 && numberOfSamplesOnHemisphere!=100 )
            {
              std::cout<<"As the current number of sampled on hemisphere is "<< numberOfSamplesOnHemisphere <<" (neither 50 nor 100), the value 50 wil be used."<<std::endl;
              numberOfSamplesOnHemisphere=50;
            }
          datFile = CreateDataFile(); // creates the file in the home directory and returns the path
        }*/
    }
//

// Added by Adrien Kaiser : test that the given nb of vertices is in the GK file
  if (numberOfSamplesOnHemisphere!=6 && numberOfSamplesOnHemisphere!=18 && numberOfSamplesOnHemisphere!=38 && numberOfSamplesOnHemisphere!=66 && numberOfSamplesOnHemisphere!=102 && numberOfSamplesOnHemisphere!=146 && numberOfSamplesOnHemisphere!=198 && numberOfSamplesOnHemisphere!=258 && numberOfSamplesOnHemisphere!=326 && numberOfSamplesOnHemisphere!=402 && numberOfSamplesOnHemisphere!=486)
    {
       std::cout<<"The number of samples on hemisphere "<<numberOfSamplesOnHemisphere<<" is not available in the file."<<std::endl;
       std::cout<<"The values available are : 6, 18 38, 66, 102, 146, 198, 258, 326, 402, 486."<<std::endl;
       std::cout<<"Abort."<<std::endl;
       exit(EXIT_FAILURE);
    }
  ODFReconstructorType ODFreconstructor ( ODFCoefsImage, numberOfSamplesOnHemisphere, numberOfSpharm, datFile );
  ODFImage = ODFreconstructor.ReconstructODFImage();

  // compute pseudo-fa image from odf
  DoubleImageType::Pointer FAImage = ODFreconstructor.GetFAImage() ;

  // write out to file if requested
  if ( writeFA ) 
    {
      std::string faFileName = output_filename_base + "_fa.nrrd" ;
      DoubleWriterType::Pointer faWriter = DoubleWriterType::New();
      faWriter->SetFileName ( faFileName );
      faWriter->SetInput( FAImage );
      faWriter->SetUseCompression(true); // Added by Adrien Kaiser

      try
	{
	  faWriter->Update();
	}
      catch ( itk::ExceptionObject & excp )
	{
	  std::cerr << "Problem writing the fa file " << faFileName << std::endl ;
	  std::cerr << excp << std::endl;
	}
    }


  // compute and write out to file pseudo - color - fa if requested
  if ( writeColorFA ) 
    {
      ColorImageType::Pointer ColorFAImage = ODFreconstructor.GetColorFAImage() ;
      std::string colorFAFileName = output_filename_base + "_color_fa.nrrd" ;
      ColorWriterType::Pointer colorFAWriter = ColorWriterType::New();
      colorFAWriter->SetFileName ( colorFAFileName );
      colorFAWriter->SetInput( ColorFAImage );
      colorFAWriter->SetUseCompression(true); // Added by Adrien Kaiser
      try
	{
	  colorFAWriter->Update();
	}
      catch ( itk::ExceptionObject & excp )
	{
	  std::cerr << "Problem writing the color fa file " << colorFAFileName << std::endl ;
	  std::cerr << excp << std::endl;
	}
    }

  if ( reconDebug ) exit ( 0 ) ;
  
  // run streamline algorithm for qc/sanity check
  if ( streamline )
    {
      ODFStreamlineType ODFStreamliner ( ODFImage, numberOfSamplesOnHemisphere * 2 ) ;
      ODFStreamliner.SetCoordinateTable ( ODFreconstructor.GetCoordinateTable () ) ;
      ODFStreamliner.SetFAImage ( FAImage ) ;
      ImageType::Pointer streamlineImage = ODFStreamliner.Streamline ( SourceImage, 100, singleTract ) ;

      std::string streamlineFileName = output_filename_base + "_streamline.nrrd" ;
      WriterType::Pointer streamlineWriter = WriterType::New();
      streamlineWriter->SetFileName ( streamlineFileName );
      streamlineWriter->SetInput( streamlineImage );
      streamlineWriter->SetUseCompression(true); // Added by Adrien Kaiser
      try
	{
	  streamlineWriter->Update();
	}
      catch ( itk::ExceptionObject & excp )
	{
	  std::cerr << "Problem writing the streamline file " << streamlineFileName << std::endl ;
	  std::cerr << excp << std::endl;
	}
    }

  // compute cost
  ODFCostType ODFCostComputer ( ODFImage, OrigFAImage, WMmaskImage, numberOfSamplesOnHemisphere * 2, alpha
				) ;
  ODFCostComputer.SetCoordinateTable ( ODFreconstructor.GetCoordinateTable () ) ;
  ODFCostComputer.SetFAImage ( FAImage ) ;
  DoubleImageType::Pointer costImage = ODFCostComputer.Cost ( SourceImage ) ;
  DoubleImageType::Pointer origImage = ODFCostComputer.GetOrigImage ( ) ;
  DoubleImageType::Pointer lengthImage = ODFCostComputer.GetLengthImage ( ) ;

  DoubleImageType::Pointer pathImage = ODFCostComputer.GetpathImage ( ) ;

  DoubleImageType::Pointer ODFVoxInImage = ODFCostComputer.GetODFVoxInImage ( ) ;
  DoubleImageType::Pointer ODFNeiOutImage = ODFCostComputer.GetODFNeiOutImage ( ) ;
  DoubleImageType::Pointer PenInOutImage = ODFCostComputer.GetPenInOutImage ( ) ;
  DoubleImageType::Pointer PenInDnImage = ODFCostComputer.GetPenInDnImage ( ) ;
  DoubleImageType::Pointer PenOutDnImage = ODFCostComputer.GetPenOutDnImage ( ) ;


  // write ODFVoxIn path image
  std::string ODFVoxInFileName = output_filename_base + "_ODFVoxIn.nrrd" ;
  DoubleWriterType::Pointer ODFVoxInWriter = DoubleWriterType::New() ;
  ODFVoxInWriter->SetFileName ( ODFVoxInFileName ); 
  ODFVoxInWriter->SetInput( ODFVoxInImage );
  ODFVoxInWriter->SetUseCompression(true); // Added by Adrien Kaiser
  try
    {
      ODFVoxInWriter->Update() ;
    }
  catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem writing the ODFVoxIn file " << ODFVoxInFileName << std::endl ;
      std::cerr << excp << std::endl ;
    }

  // write out ODFNeiOut image
  std::string ODFNeiOutFileName = output_filename_base + "_ODFNeiOut.nrrd" ;
  DoubleWriterType::Pointer ODFNeiOutWriter = DoubleWriterType::New() ;
  ODFNeiOutWriter->SetFileName ( ODFNeiOutFileName ); 
  ODFNeiOutWriter->SetInput( ODFNeiOutImage );
  ODFNeiOutWriter->SetUseCompression(true); // Added by Adrien Kaiser
  try
    {
      ODFNeiOutWriter->Update() ;
    }
  catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem writing the ODFNeiOut file " << ODFNeiOutFileName << std::endl ;
      std::cerr << excp << std::endl ;
    }

  // write out PenInOut image
  std::string PenInOutFileName = output_filename_base + "_PenInOut.nrrd" ;
  DoubleWriterType::Pointer PenInOutWriter = DoubleWriterType::New() ;
  PenInOutWriter->SetFileName ( PenInOutFileName ); 
  PenInOutWriter->SetInput( PenInOutImage );
  PenInOutWriter->SetUseCompression(true); // Added by Adrien Kaiser
  try
    {
      PenInOutWriter->Update() ;
    }
  catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem writing the PenInOut file " << PenInOutFileName << std::endl ;
      std::cerr << excp << std::endl ;
    }

  // write out PenInDn image
  std::string PenInDnFileName = output_filename_base + "_PenInDn.nrrd" ;
  DoubleWriterType::Pointer PenInDnWriter = DoubleWriterType::New() ;
  PenInDnWriter->SetFileName ( PenInDnFileName ); 
  PenInDnWriter->SetInput( PenInDnImage );
  PenInDnWriter->SetUseCompression(true); // Added by Adrien Kaiser
  try
    {
      PenInDnWriter->Update() ;
    }
  catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem writing the PenInDn file " << PenInDnFileName << std::endl ;
      std::cerr << excp << std::endl ;
    }

  // write out PenOutDn image
  std::string PenOutDnFileName = output_filename_base + "_PenOutDn.nrrd" ;
  DoubleWriterType::Pointer PenOutDnWriter = DoubleWriterType::New() ;
  PenOutDnWriter->SetFileName ( PenOutDnFileName ); 
  PenOutDnWriter->SetInput( PenOutDnImage );
  PenOutDnWriter->SetUseCompression(true); // Added by Adrien Kaiser
  try
    {
      PenOutDnWriter->Update() ;
    }
  catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem writing the PenOutDn file " << PenOutDnFileName << std::endl ;
      std::cerr << excp << std::endl ;
    }


  // write out path image
  std::string pathFileName = output_filename_base + "_path.nrrd" ;
  DoubleWriterType::Pointer pathWriter = DoubleWriterType::New() ;
  pathWriter->SetFileName ( pathFileName ); 
  pathWriter->SetInput( pathImage );
  pathWriter->SetUseCompression(true); // Added by Adrien Kaiser
  try
    {
      pathWriter->Update() ;
    }
  catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem writing the path file " << pathFileName << std::endl ;
      std::cerr << excp << std::endl ;
    }

  // write out cost image
  std::string costFileName = output_filename_base + "_cost.nrrd" ;
  DoubleWriterType::Pointer costWriter = DoubleWriterType::New() ;
  costWriter->SetFileName ( costFileName ); 
  costWriter->SetInput( costImage );
  costWriter->SetUseCompression(true); // Added by Adrien Kaiser
  try
    {
      costWriter->Update() ;
    }
  catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem writing the cost file " << costFileName << std::endl ;
      std::cerr << excp << std::endl ;
    }

  // write out orig image
  std::string origFileName = output_filename_base + "_orig.nrrd" ;
  DoubleWriterType::Pointer origWriter = DoubleWriterType::New() ;
  origWriter->SetFileName ( origFileName ); 
  origWriter->SetInput( origImage );
  origWriter->SetUseCompression(true); // Added by Adrien Kaiser
  try
    {
      origWriter->Update() ;
    }
  catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem writing the orig file " << costFileName << std::endl ;
      std::cerr << excp << std::endl ;
    }


  // write out length image
  std::string lengthFileName = output_filename_base + "_length.nrrd" ;
  DoubleWriterType::Pointer lengthWriter = DoubleWriterType::New() ;
  lengthWriter->SetFileName ( lengthFileName ); 
  lengthWriter->SetInput( lengthImage );
  lengthWriter->SetUseCompression(true); // Added by Adrien Kaiser
  try
    {
      lengthWriter->Update() ;
    }
  catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem writing the orig file " << costFileName << std::endl ;
      std::cerr << excp << std::endl ;
    }
  
  return 0;
}

