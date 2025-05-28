#include <iostream>
#include <variant>
#include "ParseArgs.h"
#include "PointSetDriver.h"
#include "FluidSPHBubbleDriver.h"
#include "LeastSquaresDriver.h"
#include "AnimationDriver.h"
#include "ArrayPointer.h"

#ifndef __Main_cpp__
#define __Main_cpp__

void Set_Threads(ParseArgs& parse_args) {
	int number_threads = parse_args.Get_Integer_Value("-tnum");
	omp_set_num_threads(number_threads);
	int max_threads = omp_get_max_threads();
	std::cout << "Set " << number_threads << " threads, run with " << max_threads << " cores\n";
}


int main(int argc,char* argv[])
{
	SparseMatrixCuda<real> mt_dev(DataHolder::DEVICE);
	SparseMatrix<real> mt_host(5, 5);
	mt_dev = mt_host;

	const int d = 3;

    ParseArgs parse_args;
    parse_args.Add_String_Argument("-o","output","output path");
    parse_args.Add_Integer_Argument("-s",400,"resolution");
	parse_args.Add_Integer_Argument("-test", 13, "test");
	parse_args.Add_Integer_Argument("-driver",2,"driver");
	parse_args.Add_Integer_Argument("-lf",10,"last frame");
	parse_args.Add_Double_Argument("-cfl",1,"cfl number");
	parse_args.Add_Integer_Argument("-tnum", 4, "Number of Threads");
	parse_args.Add_Integer_Argument("-verbose", 1, "Verbose or not");
	parse_args.Add_Integer_Argument("-diagnosis", 0, "whether to print diagnosis");
	parse_args.Add_Integer_Argument("-timing", 0, "whether output timing information");
	parse_args.Add_Integer_Argument("-from", 0, "first farme");
	parse_args.Add_String_Argument("-txt", "scene.txt", "Description File");
	parse_args.Add_Integer_Argument("-saveall", 0, "Save Mode");//if save everything
	parse_args.Add_String_Argument("-init_snapshot", "none", "initial snapshot");
    parse_args.Parse(argc,argv);

	Set_Threads(parse_args);

    std::string output_dir=parse_args.Get_String_Value("-o");
    const int scale=parse_args.Get_Integer_Value("-s");
	const int driver=parse_args.Get_Integer_Value("-driver");
	const int test=parse_args.Get_Integer_Value("-test");
	const int last_frame=parse_args.Get_Integer_Value("-lf");
	const int verbose = parse_args.Get_Integer_Value("-verbose");
	const int diagnosis = parse_args.Get_Integer_Value("-diagnosis");
	const int timing = parse_args.Get_Integer_Value("-timing");
	const int first_frame = parse_args.Get_Integer_Value("-from");
	const real cfl=parse_args.Get_Double_Value("-cfl");
	switch(driver){
	case 1:{
		PointSetDriver<d> driver;
		driver.scale=scale;
		driver.output_dir=output_dir;
		driver.test=test;
		driver.last_frame=last_frame;
		driver.cfl=cfl;
		driver.Initialize();
		driver.Run();	
	}break;
	case 2:{
		FluidSPHBubbleDriver<d> driver;
		driver.scale=scale;
		driver.output_dir=output_dir;
		driver.test=test;
		driver.last_frame=last_frame;
		driver.cfl=cfl;
		driver.verbose = verbose;
		driver.fluid.diagnosis = diagnosis;
		driver.fluid.timing = timing;
		driver.first_frame = first_frame;
		driver.scene_file_name = parse_args.Get_String_Value("-txt");
		driver.save_all = parse_args.Get_Integer_Value("-saveall");
		driver.init_snapshot_name = parse_args.Get_String_Value("-init_snapshot");
		driver.Initialize();
		driver.Run();		
		driver.fluid.output_file.close();
	}break;
	case 3:{
		LeastSquaresDriver<d> driver;
		driver.scale=scale;
		driver.output_dir=output_dir;
		driver.test=test;
		driver.last_frame=last_frame;
		driver.cfl=cfl;
		driver.Initialize();
		driver.Run();
	}break;
	case 4: {
		/**
		 * Usage:	.\point_set_surface.exe -driver 4 -test 1 -lf 30 -o wolfman
		 * Outputs:	points and meshes that can be shown in opengl_viewer
		 *			.obj files of	1) the original meshes; 
		 *							2) noisy positions;
		 *							3) voronoi graph of the original mesh 
		 */
		AnimationDriver<d> driver;
		driver.scale=scale;
		driver.output_dir=output_dir;
		driver.test=test;
		driver.last_frame=last_frame;
		driver.cfl=cfl;
		driver.Initialize();
		driver.Run();
	}break;
	}
}

#endif