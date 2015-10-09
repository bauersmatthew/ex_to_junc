/*
 * File: ex_to_junc.cpp
 * Edited: 09 Oct 2015
 * Author: Matthew Bauer
 */

#include <string>
#include <vector>
#include <stdint.h>
#include <fstream>
#include <sstream>
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "api/SamHeader.h"
#include <cstdio>

using namespace BamTools;

const int ARGC_MIN = 3;
const int ARGC_MAX = 3;

std::string Eify(const std::string& msg);
void PrintUsage();

struct Sec
{
	uint32_t start;
	uint32_t end;
};
std::vector<Sec> ReadBed(const std::string& path);
Sec GetSecFromBedLine(const std::string& line);
std::vector<std::string> SepTab(const std::string& str);
template<typename T>
bool StrToT(const std::string& str, T *ptr);

const std::string TMP_PATH = "ex_to_junc_tmp_bamout_file_28lkali2luao8roiuwer-ewrwehrlwr0";
void RMPBam(BamReader&, std::vector<Sec>&);

int main(int argc, char **argv)
{
	// check for --help param
	if(argc == 2 && std::string(argv[1]) == "--help")
	{
		PrintUsage();
		return 0;
	}
	
	// check usage otherwise
	if(argc < ARGC_MIN || argc > ARGC_MAX)
	{
		std::cerr << Eify("incorrect usage") << std::endl;
		PrintUsage();
		return 1;
	}
	
	// read bed
	std::vector<Sec> secs;
	try
	{
		secs = ReadBed(argv[2]);
	}
	catch(std::string& se)
	{
		std::cerr << Eify(se) << std::endl;
		return 2;
	}
	
	// read/process/write bam
	BamReader br;
	if(!br.Open(argv[1]))
	{
		std::cerr << Eify("failed to open file '" + std::string(argv[1]) + "'") << std::endl;
		return 3;
	}
	try
	{
		RMPBam(br, secs);
	}
	catch(std::string& se)
	{
		br.Close();
		std::cerr << Eify(se) << std::endl;
		return 4;
	}
	br.Close();
	
	// we should be done
	return 0;
}

std::string Eify(const std::string& msg)
{
	return "e: " + msg;
}
void PrintUsage()
{
	std::cerr << "Usage: ex_to_junc <bam> <circ-bed>" << std::endl;
}

std::vector<Sec> ReadBed(const std::string& path)
{
	// open file
	std::ifstream fin(path);
	if(!fin)
	{
		throw "failed to open file '" + path + "'";
	}
	
	std::vector<Sec> ret;
	Sec csec;

	// read line by line
	std::string line;
	std::getline(fin, line);
	while(!fin.eof())
	{
		if(line.empty())
		{
			// skip
			std::getline(fin, line);
			continue;
		}
		
		// convert to sec
		try
		{
			csec = GetSecFromBedLine(line);
		}
		catch(...)
		{
			fin.close();
			throw;
		}
		
		// add to vector, get next line
		ret.push_back(csec);
		std::getline(fin, line);
	}
	fin.close();
	
	return ret;
}
Sec GetSecFromBedLine(const std::string& line)
{
	// break up
	std::vector<std::string> seps = SepTab(line);
	if(seps.size() < 3) // bed format: chrom, start, end, [name]
	{
		throw "invalid bed line";
	}
	
	// convert
	Sec ret;
	if(!StrToT<uint32_t>(seps.at(1), &ret.start))
		throw "invalid bed line";
	if(!StrToT<uint32_t>(seps.at(2), &ret.end))
		throw "invalid bed line";
	
	// return
	return ret;
}
std::vector<std::string> SepTab(const std::string& str)
{
	std::string cs = "";
	std::vector<std::string> ret;
	for(int i = 0; i < str.size(); i++)
	{
		char ch = str.at(i);
		if(ch == '\t')
		{
			ret.push_back(cs);
			cs = "";
		}
		else
		{
			cs += ch;
		}
	}
	ret.push_back(cs);
	return ret;
}
template<typename T>
bool StrToT(const std::string& str, T *ptr)
{
	std::stringstream ss;
	T t;
	ss << str;
	if(!(ss >> t))
		return false;
	if(ptr != NULL)
		*ptr = t;
	return true;
}

void RMPBam(BamReader& br, std::vector<Sec>& secs)
{
	// check for existance of tmp file
	{
		std::ifstream tmpf_check(TMP_PATH);
		if(tmpf_check)
		{
			tmpf_check.close();
			throw "temp file already exists -- please rename or delete";
		}
	}

	// open writer
	BamWriter bw;
	if(!bw.Open(TMP_PATH, br.GetHeader(), br.GetReferenceData()))
	{
		throw "could not open temp file for writing";
	}

	BamAlignment baln;
	while(true)
	{
		if(!br.GetNextAlignment(baln))
		{
			// probably means that we have read all of them
			break;
		}
		if(!baln.IsMapped())
		{
			// regurgitate
			bw.SaveAlignment(baln);
			continue;
		}
		
		// search for it in the bed database
		uint32_t start_pos = baln.Position;
		uint32_t end_pos = baln.Position + baln.Length;
		bool found = false;
		for(uint64_t i = 0; i < secs.size(); i++)
		{
			Sec csec = secs.at(i);
			
			// check for mapped to the actual circ region at all
			if(start_pos < csec.start || end_pos > csec.end)
			{
				// nope
				continue;
			}
			
			// check for mapped across the junction
			uint32_t csec_midpos = (csec.end - csec.start)/2 + csec.start;
			if(start_pos > csec_midpos || end_pos < csec_midpos)
			{
				// nope. but this time there's no chance of it working anywhere else.
				break;
			}
			
			// all the reqs have been met; it is a valid alignment.
			found = true;
			break;
		}
		
		// process it being found or not
		if(found)
		{
			// regurgitate
			bw.SaveAlignment(baln);
		}
		else
		{
			// set as not mapped, then write
			baln.SetIsMapped(false);
			bw.SaveAlignment(baln);
		}
	}
	bw.Close(); // we are done writing now. yay.
	
	// read it back in, output to the console while doing so
	std::ifstream fin(TMP_PATH, std::ios::in | std::ios::binary);
	if(!fin)
	{
		// I don't know what the fuck would cause this
		throw "could not read temp file after writing it";
	}
	char read_char = fin.get();
	while(!fin.eof())
	{
		std::cout << read_char;
		read_char = fin.get();
	}
	fin.close();
	
	// delete the temp file
	remove(TMP_PATH.c_str());
	
	// done!
}
