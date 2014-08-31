// LiquidBlender.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "gtest/gtest.h"
#include "ICP.hpp"
#include "Animation.h"
#include "Utility.hpp"
#include "Blend.hpp"
#include "Configuration.hpp"
#include "Output.hpp"
#include <Bench/BenchTimer.h>
#include "Retiming.h"
#include "Retimer.hpp"

Float Globals::AvgEdgeLength = 0;
const int Globals::InvalidNeighbour = -1;	

Eigen::IOFormat CommaFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "  ", ";");

void setupGoogleLogging() {
	
	for (int i = 0; i < 4; i++)
		google::SetLogDestination(i, "logs/log");
}


template<int DIM>
void runICP(Configuration<DIM>& cfg) {
	ICP<DIM> icp(cfg);

	icp.initialize();
			
	for (int i = 1; i <= cfg.maxIterations; i++) {		
		icp.doIteration();			
	}

	icp.cleanup();

	blend(cfg, icp.getCorrAnim(), 50);
}

template<int DIM>
void retimeAnimation(RetimingConfiguration<DIM>& cfg) {

}

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
	char ** itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
	{
		return *itr;
	}
	return 0;
}



bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
	return std::find(begin, end, option) != end;
}

int main(int argc, char* argv[])
{	
	std::string progName(argv[0]);
	std::string str( progName.begin(), progName.end() );
	google::InitGoogleLogging(str.c_str());
	setupGoogleLogging();

	Eigen::BenchTimer timer;
	timer.reset();
	timer.start();

	if (cmdOptionExists(argv, argv+argc, "--test")) {
		::testing::InitGoogleTest(&argc, argv);
		return RUN_ALL_TESTS();
	}

	int dim = -1;

	if(cmdOptionExists(argv, argv+argc, "--dim"))
	{
		std::string dimString = getCmdOption(argv,argv+argc, "--dim");
		dim = stoi(dimString);
	}
	else {
		LOG(ERROR) << "Dimension needs to be specified using --dim" << std::endl;
		return -1;
	}

	if(cmdOptionExists(argv, argv+argc, "--register"))
	{
		std::string wCfgFile = getCmdOption(argv,argv+argc, "--register");
		std::string cfgFile( wCfgFile.begin(), wCfgFile.end() );			

		if (dim == 2) {
			Configuration<2> cfg;
			if (!cfg.load(cfgFile)) {
				LOG(ERROR) << "Failed to load configuration file:" << cfgFile << std::endl;
				return -1;
			}

			runICP<2>(cfg);
		}
		else if (dim == 3) {
			Configuration<3> cfg;
			if (!cfg.load(cfgFile)) {
				LOG(ERROR) << "Failed to load configuration file:" << cfgFile << std::endl;
				return -1;
			}

			runICP<3>(cfg);
		}				
	}

	if (cmdOptionExists(argv, argv+argc, "--blend")) {
		std::string corrFile;
		std::string blendString;
		int blendValue;


		if(cmdOptionExists(argv, argv+argc, "--corr")) {
			 corrFile = getCmdOption(argv,argv+argc, "--corr");
		}
		else {
			LOG(ERROR)<<"Please specify correspondence file using --corr"<<std::endl;
			return -1;
		}
		if(cmdOptionExists(argv, argv+argc, "--blendValue")) {
			blendString = getCmdOption(argv,argv+argc, "--blendValue");
			blendValue = std::stoi(blendString);
		}
		else {
			LOG(ERROR)<<"Please blend percentange using --blendValue"<<std::endl;
			return -1;
		}
				
		LOG(INFO) << "Blending :" << corrFile.c_str() << " " << blendValue << std::endl;
		
		if ((dim == 2 && !blend<2> ( corrFile, blendValue)) ||
			(dim == 3 && !blend<3> ( corrFile, blendValue))
			)
		{
			LOG(ERROR) << "Failed to blend animation!" << std::endl;
			return -1;
		}
	}
	else if (cmdOptionExists(argv, argv+argc, "--blendBary")) {

		string corrFiles[2];
		Float baryU, baryV, baryW;

		if(cmdOptionExists(argv, argv+argc, "--corr0")) {
			corrFiles[0] = getCmdOption(argv,argv+argc, "--corr0");
		}		
		if(cmdOptionExists(argv, argv+argc, "--corr1")) {
			corrFiles[1] = getCmdOption(argv,argv+argc, "--corr1");
		}		

		if(cmdOptionExists(argv, argv+argc, "--U")) {
			string tmp;
			tmp = getCmdOption(argv,argv+argc, "--U");
			baryU = std::atof(tmp.c_str());
		}		
		if(cmdOptionExists(argv, argv+argc, "--V")) {
			string tmp;
			tmp = getCmdOption(argv,argv+argc, "--V");
			baryV = std::atof(tmp.c_str());
		}		
		if(cmdOptionExists(argv, argv+argc, "--W")) {
			string tmp;
			tmp = getCmdOption(argv,argv+argc, "--W");
			baryW = std::atof(tmp.c_str());
		}		

		if (dim == 3 && !blend<3> ( corrFiles, baryU, baryV, baryW) )
		{
			LOG(ERROR) << "Failed to blend animation!" << std::endl;
			return -1;
		}
	}
	else if (cmdOptionExists(argv, argv+argc, "--retime")) {
		std::string wCfgFile = getCmdOption(argv,argv+argc, "--retime");
		std::string cfgFile( wCfgFile.begin(), wCfgFile.end() );			

		if (dim == 2) {
			RetimingConfiguration<2> cfg;
			if (!cfg.load(cfgFile)) {
				LOG(ERROR)  << "Failed to load configuration file:" << cfgFile << std::endl;
				return -1;
			}

			Retimer<2> retimer(cfg);
			
		}
		else if (dim == 3) {
			RetimingConfiguration<3> cfg;
			if (!cfg.load(cfgFile)) {
				LOG(ERROR)  << "Failed to load configuration file:" << cfgFile << std::endl;
				return -1;
			}

			Retimer<3> retimer(cfg);

		}
	}


	
	timer.stop();
	LOG(INFO) << "Execution time: " << timer.total() << endl;
}

