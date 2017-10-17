/*
	vAlign is a fast verifier of reads of aligned DNA sequence,
	recalled from initial artificial FastQ sequence.
	It compares the original and actual coordinates of each read
	and prints statistics of right and wrong mappings.
	
	Copyright (C) 2017 Fedor Naumenko (fedor.naumenko@gmail.com)

	This program is free software. It is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY;
	without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
	See the	GNU General Public License for more details.
 */

#include "Data.h"
#include "vAlign.h"
#include <fstream>

using namespace std;

const string Product::Title = "vAlign";
const string Product::Version = "1.0";
const string Product::Descr = "verify Alignment";

const string OutFile = string(Product::Title) +  "_out.txt";
const string HelpOutFile = "duplicate standard output to " + OutFile + " file";

enum eOptGroup	{ oINPUT, oOUTPUT, oOTHER };	// oOTHER should be the last 
const BYTE	Options::_GroupCount = oOTHER + 1;	// count of option groups in help

const char* Options::_OptGroups [] = {			// names of option groups in help
	"Input", "Output", "Other"
};

// --info option: types of info notations
const char* infos [] = { "NM", "CNT", "STAT" };	// corresponds to eInfo; iNONE and iLAC are hidden

//	{ char,	str,	Signs,	type,	group,	defVal,	minVal,	maxVal,	strVal,	descr, addDescr }
// field 7: vUNDEF if value is prohibited
// field 6: vUNDEF if no default value should be printed
Options::Option Options::_Options [] = {
	{ 'g', "gen",	1, tNAME, oINPUT, vUNDEF, 0, 0, NULL,
	"reference genome library or single nucleotide sequence.", NULL },
	{ 'c', Chrom::Abbr,	0,	tCHAR,	oINPUT, vUNDEF, 0, 0, NULL,	"treat specified chromosome only", NULL },
	{ HPH, "min-scr",	0, tINT,	oINPUT, vUNDEF, 0, 1000, NULL, "score threshold for treated reads", NULL },
	{ HPH, "char-case",	0,	tENUM,	oINPUT, FALSE,	0, 2, (char*)Options::Booleans,
	"recognize uppercase and lowercase characters in template and test\nas different", NULL },
	{ 'i', "info",	0,	tENUM, oOUTPUT,	Obj::iSTAT, Obj::iNM, Obj::iSTAT, (char*)infos,
	"print information about file:\n? - name only, ? - number of reads, ? - statistics", NULL },
	{ 'w', "warn",	0,	tENUM,	oOUTPUT, FALSE,	vUNDEF, 2, NULL,
	"print each read ambiguity, if they exist" },
	{ 'o', "out",	0,	tENUM,	oOUTPUT,FALSE,	vUNDEF, 2, NULL, HelpOutFile.c_str(), NULL },
	{ 't', "time",	0,	tENUM,	oOTHER,	FALSE,	vUNDEF, 2, NULL, "print run time", NULL },
	{ 'v', Version,	0,	tVERS,	oOTHER,	vUNDEF, vUNDEF, 0, NULL, "print program's version", NULL },
	{ 'h', "help",	0,	tHELP,	oOTHER,	vUNDEF, vUNDEF, 0, NULL, "print usage information", NULL }
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
	if (argc < 2)	return Options::PrintUsage(false);			// output tip
	int fileInd = Options::Tokenize(argc, argv);
	if( fileInd < 0 )	return 1;								// wrong option
	if(!Chrom::SetStatedID(Options::GetSVal(oCHROM))) return 1;	// wrong chrom name

	int ret = 0;						// main() return code
	if( Options::GetBVal(oOUTFILE) )	outfile.open(OutFile.c_str());
	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		// check file names first of all
		const char* gName = FS::CheckedFileDirName(oGFILE);
		const char* aName = FS::CheckedFileName(argv[fileInd]);
		ChromFiles cFiles(FS::CheckedFileDirName(oGFILE));
		ChromSizes cSizes(cFiles);

		BedR test("test alignment", aName, &cSizes, 
			Obj::eInfo (Options::GetIVal(oINFO)), 
			true, true, Options::GetBVal(oALARM),
			true,		Options::GetIVal(oMINSCR));
		if( test.ReadNameType() != Read::nmPos )
			Err("Read's name does not contain read's initial position").Throw();

		vAlign(cFiles, test);
	}
	catch(Err &e)		{ ret = 1;	dout << e.what() << EOL; }
	catch(exception &e)	{ ret = 1;	dout << e.what() << EOL; }
	catch(...)			{ ret = 1;	dout << "Unregistered error\n"; }

	timer.Stop();
	if( outfile.is_open() )		outfile.close();	// in case of holding execution by user
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
		if(!Chrom::StatedAll && Chrom::StatedID() != cID)	continue;
		dout << Chrom::AbbrName(cID) << EOL;
		_preciseAccum.Reset();
		_mismAccums.Clear();
		Nts nts(cFiles.FileName(cID), false);
		ritend = bedR.ReadsEnd(cit);
		for(rit=bedR.ReadsBegin(cit); rit!=ritend; rit++)
			if( rit->InitCID == cID )
				if( rit->Pos==rit->Num )	// Pos: actual start position, Num: initial start position
					_preciseAccum.AddRead(rit->Score);
				else
					_mismAccums[ VerifyRead(nts, rit->Num, rit->Pos) ].AddRead(rit->Score);
		PrintStats(cID, bedR.ReadsCount(cID));
	}
}

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
	dout << "total reads per chrom " << Chrom::Name(cID) << ":\t" << cnt << TAB
		 << sPercent(cnt, rCnt, 0, 0, false) //<< EOL;
		 << TAB << rCnt << EOL;
	cnt = rCnt - cnt;
	dout << "reads per different chroms:\t" << cnt << TAB << sPercent(cnt, rCnt, 0, 0, false) << EOL;
}
/************************ end of class vAlign ************************/
