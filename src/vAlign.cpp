/*
 * verifying alignment
 * Copyright (C) 2016 Fedor Naumenko
 */

#include "Data.h"
#include "vAlign.h"
#include <fstream>

using namespace std;

const string Product::Title = "vAlign";
const string Product::Version = "1.0";
const string Product::Descr = "verify Alignment";

const string Tested = "tested";
const string Alignment = " alignment";
const string OutFile = string(Product::Title) +  "_out.txt";
const string HelpOutFile = "duplicate standard output to " + OutFile + " file";

enum eOptGroup	{ oINPUT, oOUTPUT, oOTHER };	// oOTHER should be the last 
const BYTE	Options::_GroupCount = oOTHER + 1;	// count of option groups in help

const char* Options::_OptGroups [] = {			// names of option groups in help
	"Input", "Ambig output", "Other"
};

Options::Option Options::_Options [] = {
	{ 'g', "gen",	1, true, tNAME, oINPUT, vUNDEF, 0, 0, NULL,
	"genome size file, or genome library, or single nucleotide sequence" },
	//{ 'T', "templ",	0, true, tNAME, oINPUT, vUNDEF, 0, 0, NULL,
	//"template sequence" },
	{ 'c',Chrom::Abbr,0,true,tCHAR, oINPUT, vUNDEF, 0, 0, NULL,
	"treat stated chromosome only (all)" },
	{ HPH, "min-scr",	0, true, tINT, oINPUT, vUNDEF, 0, 1000, NULL, "score threshold for treated reads (lack)" },
	{ HPH, "char-case",0,true, tENUM,	oINPUT, FALSE,	0, 2, (char*)Options::Booleans,
	"recognize uppercase and lowercase characters in template and test\nas different" },
	{ HPH, "alarm",	0, false,tENUM, oOUTPUT,FALSE,	0, 2, NULL,
	"output features ambiguities, if they exist" },
	{ HPH, "stat",	0, false,tENUM, oOUTPUT,FALSE,	0, 2, NULL,
	"output features ambiguities statistics, if they exist" },
	{ 'o', "out",	0, false,tENUM, oOUTPUT,FALSE,	0, 2, NULL, HelpOutFile.c_str() },
	{ 't', "time",	0, false,tENUM, oOTHER,	FALSE,	0, 2, NULL, "output run time" },
	{ 'v', Version,	0, false,tVERS,oOTHER, vUNDEF, 0, 0, NULL, "print program's version and quit" },
	{ 'h', "help",	0, false,tHELP,oOTHER, vUNDEF, 0, 0, NULL, "print usage information and quit" }
};

const BYTE	Options::_OptCount = oHELP + 1;
const BYTE	Options::_UsageCount = 1;		// count of 'Usage' variants in help
const Options::Usage Options::_Usages[] = {	// content of 'Usage' variants in help
	{	vUNDEF,	" sequence"	}
};

ofstream outfile;		// file ostream duplicated cout; inizialised by file in code
dostream dout(cout, outfile);	// stream's duplicator

/*****************************************/
int main(int argc, char* argv[])
{
	//const char* a = "12345-6789/2";
	//long b = atol(a+6);
	//cout << b << EOL;	return 0;
	if (argc < 2)	
	{ Options::PrintUsage(false);	return 0; }		// output tip

		int fileInd = Options::Tokenize(argc, argv, (Tested + Alignment).c_str());
	if( fileInd < 0 )	return 1;		// wrong option
	int ret = 0;						// main() return code

	bool getTime = Options::GetBVal(oTIME);
	if( Options::GetBVal(oOUTFILE) )	outfile.open(OutFile.c_str());
	Timer timer(getTime);
	timer.Start();
	try {
		// check file names first of all
		const char* genName	  = FS::CheckedFileDirName(oGFILE);
		//const char* templName = FS::CheckedFileName(oTEMPL);
		const char* testlName = FS::CheckedFileName(argv[fileInd]);
		bool alarm	= Options::GetBVal(oALARM);		// print ambiguity messages
		bool stats	= Options::GetBVal(oSTATS);		// print ambiguity statistics
		chrid cID	= Chrom::ID(Options::GetSVal(oCHROM));
		
		BedR test((Tested + Alignment + MSGSEP_BLANK).c_str(),
			testlName, cID, true, getTime, true, alarm, stats, true, true, Options::GetIVal(oMINSCR));
		if( test.ReadNameType() != Read::nmPos )
			Err("Read's name does not contain read's initial position").Throw();

		Timer timer(getTime);
		timer.Start();
		vAlign(ChromFiles(genName, cID), test);
		timer.Stop();
	}
	catch(Err &e)		{ ret = 1;	dout << e.what() << EOL; }
	catch(exception &e)	{ ret = 1;	dout << e.what() << EOL; }
	catch(...)			{ ret = 1;	dout << "Unregistered error\n"; }

	timer.Stop("total: ", false, true);
	if( outfile.is_open() )		outfile.close();	// in case of holding execution by user
//#ifdef OS_Windows
//	system("pause");
//#endif
	return ret;
}

/************************ class vAlign ************************/

vAlign::vAlign(const ChromFiles& cFiles, BedR &bedR) :
	_caseDiff(Options::GetBVal(oCCASE)),
	_maxScore(bedR.MaxScore())
{
	BedR::cItemsIter rit, ritend;
	chrid	cID;

	Read::Len = bedR.ReadLen();
	_mismAccums.Reserve(Read::Len+1);
	for(BedR::cIter cit=bedR.cBegin(); cit!=bedR.cEnd(); cit++) {
		cID = CID(cit);
		dout << Chrom::AbbrName(cID) << EOL;
		_preciseAccum.Reset();
		_mismAccums.Clear();
		Nts nts(cFiles.FileName(cID), false);
		ritend = bedR.ReadsEnd(cit);
		for(rit=bedR.ReadsBegin(cit); rit!=ritend; rit++)
			if( rit->InitCID == cID )
				if( rit->Pos==rit->Num )
					_preciseAccum.AddRead(rit->Score);
				else
					_mismAccums[ VerifyRead(nts, rit->Num, rit->Pos) ].AddRead(rit->Score);
		PrintStats(cID, bedR.ReadsCount(cID));
	}
}

//vAlign::vAlign(const ChromFiles& cFiles, BedR &bedRtempl, BedR &bedRTest)
//{
//	_caseDiff = Options::GetBVal(oCCASE);
//	_maxScore = bedRTest.MaxScore();
//	BedR::ItemsIter rit1, rit2, itr1begin, itr2begin, rit1end, rit2end;
//	BedR::cIter cit2;
//	readlen errCnt;
//
//	if( (Read::Len = bedRtempl.ReadLen()) != bedRTest.ReadLen() )
//		Err("different length of read").Throw();
//	_mismAccums.Init(Read::Len+1);
//
//	Timer timer(true);
//	timer.Start();
//	for(BedR::cIter cit1=bedRtempl.cBegin(); cit1!=bedRtempl.cEnd(); cit1++) {
//		if( !bedRTest.FindChrom(CID(cit1)) )	continue;
//		
//		dout << Chrom::AbbrName(CID(cit1)) << EOL;
//		Nts nts( cFiles.FileName(CID(cit1)), 0, true, false);
//		cit2 = bedRTest.GetIter(CID(cit1));
//		bedRtempl.SortByNumb(cit1);
//		bedRTest.SortByNumb(cit2);
//		
//		itr1begin = bedRtempl.ReadsBegin(cit2);
//		itr2begin = bedRTest.ReadsBegin(cit2);
//		rit1end = bedRtempl.ItemsEnd(cit1);
//		rit2end = bedRTest.ItemsEnd(cit2);
//
//		for(rit1=itr1begin; rit1!=rit1end; rit1++)
//			for(rit2=itr2begin; rit2!=rit2end; rit2++) {
//				if( rit1->Numb < rit2->Numb)	break;
//				if( rit1->Numb == rit2->Numb ) {
//					errCnt = rit1->Pos==rit2->Pos ? 0 : VerifyRead(nts, rit1->Pos, rit2->Pos);
//					_mismAccums[errCnt].AddRead(rit2->Score);
//					itr2begin = rit2;
//					break;
//				}
//			}
//	}
//}

// Gets count of mismatches for tested Read
//	@nts: chromosome sequence
//	@templPos: start position of template Read
//	@testPos: start position of tested Read
//	return: count of testet Read's mismatches in comparison with template Read
readlen vAlign::VerifyRead(const Nts& nts, chrlen templPos, chrlen testPos)
{
	const char* pos1 = nts.Read(templPos);
	const char* pos2 = nts.Read(testPos);
	readlen cnt = 0;

	for(readlen i=0; i<Read::Len; i++) {
		if( _caseDiff ) {
			if(*(pos1+i) != *(pos2+i))
				cnt++;
		}
		else
			if(toupper(*(pos1+i)) != toupper(*(pos2+i)))
				cnt++;
	}
	return cnt;
}

// Prints statistic for given chrom
//	@cID: chrom ID
//	@rCnt: total count of Reads for given chrom
void vAlign::PrintStats(chrid cID, size_t rCnt)
{
	dout << "mism\treadCnt\tquality\n";
	dout << "precise\t";	_preciseAccum.Print(_maxScore);
	//dout << "_lowScoreCnt = " << _lowScoreCnt << EOL;
	size_t cnt = _preciseAccum.Count();		// count of Reads at given chrom
	for(readlen i=0; i<=Read::Len; i++) {
		dout << int(i) << TAB;
		_mismAccums[i].Print(_maxScore);
		cnt += _mismAccums[i].Count();
	}
	dout << "total reads at chrom " << Chrom::Name(cID) << ":\t" << cnt << TAB
		 << sPercent(cnt, rCnt, 0, 0, false) //<< EOL;
		 << TAB << rCnt << EOL;
	cnt = rCnt - cnt;
	dout << "reads at different chrom:\t" << cnt << TAB << sPercent(cnt, rCnt, 0, 0, false) << EOL;
}
/************************ end of class vAlign ************************/
