#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <math.h>
#include "MyMath.h"
#include "ExpData.h"
#include <stdlib.h>

using namespace std;

ExpData::ExpData()
{
	mSize=-1;
}
ExpData::~ExpData()
{
}

bool ExpData::getExpDataFromFile(char *filename)
{
	ifstream infile(filename);
	if( !infile.is_open() ) {
		cerr << "cannot open file: " << filename << endl;
		return false;
	}

	string line;
	int nline = 0;
//	double max_value = -1.0;
	while( !infile.eof() ) {
		getline( infile, line );
		if( line.length() == 0 ) continue;
		if( line[line.length()-1] == '\r' ) { line.erase(line.length()-1) = '\0'; }
		nline++;
		if( line[0] == '#' ) continue;

		// fast reading
		char str[10240];
		//assert(line.length() < sizeof(str));
		strcpy(str,line.c_str());
		char *p[1024];
		int np = 1, n=strlen(str);
		p[0] = str;
		for( int i=0 ; i<n ; i++ ) {
			if( str[i] == '\t' ) { str[i] = '\0'; p[np] = &str[i+1]; np++; }
		}
		int ntokens = np;
		mSize=ntokens;
		// check the number of items

/*
		int d = ntokens/2;
		if( ntokens%2 == 1 ) d--;
		if( mSize == 0 || d != mSize ) {
			cerr << "input file error: not enough fileds" << endl;
			cerr << filename << "(" << nline << "): " << line << endl;
			infile.close();
			return false;
		}
*/
		// fill data
		expdata_t expdata;
		expdata.name = p[0];
		for( int i=1 ; i<ntokens ; i++ ) {

			expdata.exp.push_back(ROUND(atof(p[i])));

			/*
			cendata.x = ROUND(atof(p[2*i]));
			cendata.e = atoi(p[2*i+1]);
			if( cendata.e != 0 ) cendata.e = 1;
			mData[i].push_back(cendata);
			if( cendata.x > max_value ) max_value = cendata.x;
			*/
		}
		mEXP.push_back(expdata);


	}

	infile.close();

	// set psuedo infinity
	//pINF = round( 10.0 * max_value );

	// cout << "read " << nline << " samples from the file in " << mDimension << " dimension." << endl;

	return true;
}
