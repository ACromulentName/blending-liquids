#include "stdafx.h"
#include "json.h"
#include "Utility.hpp"
#include "gtest/gtest.h"

TEST(JSONTest, ReadFile) {
	Json::Value root;   // will contains the root value after parsing.
	Json::Reader reader;	
	std::string fileContents;
	std::string filename = "configs/configTemplate.json";
	ASSERT_TRUE(readFileIntoString(filename, fileContents));

	bool parsingSuccessful = reader.parse(fileContents, root );

	EXPECT_TRUE(parsingSuccessful);
	if ( !parsingSuccessful )
	{
		// report to the user the failure and their locations in the document.
		std::cout  << "Failed to parse configuration\n"
			<< reader.getFormattedErrorMessages();
		return;
	}

	const Json::Value animA = root["animA"];
	if (!animA.isNull()) {
		EXPECT_TRUE(animA["path"].asString()=="pathA");
	}
}

TEST(JSONTest, ReadConfig) {
	Json::Value root;   // will contains the root value after parsing.
	Json::Reader reader;	
	std::string fileContents;
	std::string filename = "configs/configTemplate.json";
	ASSERT_TRUE(readFileIntoString(filename, fileContents));

	bool parsingSuccessful = reader.parse(fileContents, root );

	EXPECT_TRUE(parsingSuccessful);
	if ( !parsingSuccessful )
	{
		// report to the user the failure and their locations in the document.
		std::cout  << "Failed to parse configuration\n"
			<< reader.getFormattedErrorMessages();
		return;
	}
	
	if (root["animA"].isNull()) {
		LOG(ERROR) << "No entry for animation A!";
		return;
	}
	else
	{
		std::string path = root["animA"]["path"].asString();
		int startFrame = root["animA"]["startFrame"].asInt();
		int endFrame = root["animA"]["endFrame"].asInt();
	}
	
	if (root["animB"].isNull()) {
		LOG(ERROR) << "No entry for animation B!";
		return;
	}
	else
	{
		std::string path = root["animB"]["path"].asString();
		int startFrame = root["animB"]["startFrame"].asInt();
		int endFrame = root["animB"]["endFrame"].asInt();
	}

	int numIterations = root.get("numIterations", 10).asInt();

	if (root["weights"].isNull()) {
		LOG(ERROR) << "No weights specified!";
	}
	else {
		Float weightRigid = root["weights"]["rigid"].asFloat();
		Float weightSmooth = root["weights"]["smooth"].asFloat();
		Float weightPoint = root["weights"]["point"].asFloat();
		Float weightPlane = root["weights"]["plane"].asFloat();
	}

	if (root["correspondences"].isNull()) {
		LOG(ERROR) << "Correspondence file is missing!";
	}
	else {
		std::string path = root["correspondences"]["path"].asString();
		bool isEnabled = root["correspondences"]["enabled"].asBool();
	}
}

TEST(JSONTest, ReadUserCorrs) {
	Json::Value root;   // will contains the root value after parsing.
	Json::Reader reader;	
	std::string fileContents;
	std::string filename = "configs/ucorrTemplate.json";
	ASSERT_TRUE(readFileIntoString(filename, fileContents));

	bool parsingSuccessful = reader.parse(fileContents, root );

	EXPECT_TRUE(parsingSuccessful);
	if ( !parsingSuccessful )
	{
		// report to the user the failure and their locations in the document.
		std::cout  << "Failed to parse configuration\n"
			<< reader.getFormattedErrorMessages();
		return;
	}

	Json::Value linesA = root["linesA"];
	if (linesA.isNull()) {
		LOG(ERROR) << "No line features in animation A!";
	}
	else {
		int numLines = linesA.size();
		for (int i = 0; i < numLines; i++) {
			std:: cout << "Line:" << i << " : ";
			for (int j = 0; j < linesA[i]["indices"].size(); j++) {				
				std::cout << linesA[i]["indices"][j].asInt() << " ";
			}
			std::cout << std::endl;
		}
	}
}