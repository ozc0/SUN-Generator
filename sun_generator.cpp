#include <htslib/faidx.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <bits/stdc++.h>

using namespace std;

#define ft first
#define sc second
#define pii pair<int,int>
#define vs vector <string>
#define pss pair<segment,segment>
#define mp make_pair
#define pb push_back

struct files_str {
  char *fasta_file;
  char *duplication_file;
  char *output_file;
} files;

struct segment {
  string chr;
  int beg, en, strand;
  segment ( string chr0, int beg0, int en0, int strand0 ) :chr(chr0), beg(beg0), en(en0), strand(strand0){}
  segment() {chr = ""; beg = en = strand = 0;}
  bool operator== ( const segment& r) const {
    if ( this->chr.compare(r.chr) == 0 && this->beg == r.beg && this->en == r.en && this->strand == r.strand )
      return true;
    return false;
  }
  bool operator<(  const segment& r ) const {
    if ( this->chr.compare(r.chr) < 0 )
      return true;
    if ( this->chr.compare(r.chr) > 0 )
      return false;
    if ( this->beg != r.beg )
      return this->beg < r.beg;
    if ( this->en != r.en )
      return this->en < r.en;
    return this->strand < r.strand;
  }
  segment& operator = (const segment& other) {
    this->chr = other.chr;
    this->beg = other.beg;
    this->en = other.en;
    this->strand = other.strand;
    return *this;
  }
  
};

struct similar {
	segment seg1, seg2;
	int after_first, len;
};

set < segment > all_segments, super_segments, suns;
map < segment, vector <segment> > close_segments;
map < segment, set < segment > > group_members;
map < segment, segment > grp;
vector < pss > pairs;
vector < similar > similars;

faidx_t* ref_fai;
FILE *save;
int uneven_cnt;

void print_help () {
  fprintf(stderr,"To run this program, you need fasta file and duplication file. Also, the output file must be specified in the following format.\n");
  fprintf(stderr,"\t-f [FASTA file]        : Reference genome in FASTA format\n");
  fprintf(stderr,"\t-d [duplication file]  : Duplication file file in tab format. (see example file for identation)\n");
  fprintf(stderr,"\t-o [output file]       : output file\n");
  fprintf(stderr,"\t-h                     : help option\n");
  fprintf(stderr,"\nExample usage: ./sun_gen -f example.fasta -d example_duplication.tab -o example.out\n");
}

void get_file_name ( int argc , char** argv ) {

  int opt;
  bool flag_f = 0, flag_d = 0, flag_o = 0;
  
  while ( ( opt = getopt( argc, argv , "f:d:o:h:" ) ) != -1 ) {
    switch ( opt ) {
    case 'f':
      files.fasta_file = optarg;
      flag_f = 1;
      break;
    case 'd':
      files.duplication_file = optarg;
      flag_d = 1;
      break;
    case 'o':
      files.output_file = optarg;
      flag_o = 1;
      break;
    case 'h':
      print_help();
      exit(0);
      break;
    default:
      fprintf(stderr, "Invalid input type\n");
      print_help();
      exit(EXIT_FAILURE);
      break;
    }
  }

  if ( !( flag_f && flag_d && flag_o ) ) {
    fprintf(stderr,"Missing one of the input/output files.\n");
    print_help();
    exit(0);
  }
}

vector < pss > read_tab();

void compare ( segment seg1, segment seg2, int aft, int len ) {

  char *ref_seq;

  // for ( int i = 0 ; i < faidx_nseq( ref_fai ) ; i++ ) {
    
  //   const char *name = faidx_iseq(ref_fai,i);


  //   int sqlen = faidx_seq_len ( ref_fai , name );
  //   int reallen;
  //   char *seq = fai_fetch ( ref_fai , name , &reallen )  ;

  //   if ( name[3] == 'X' || name[3] == 'Y' || name[3] == 'M' ) {
  //     printf(">%s\n",name);
  //     for ( int j = 0 ; j < reallen ; j++ ) {
		// 		if (  j%50 == 0 && j != 0 )
	 //  			puts("");
		// 		printf("%c",seq[j]);
  //     }
  //     puts("");
  //   }
    
  //   else {
  //     printf(">%sp1\n",name);
  //     for ( int j = 0 ; j < reallen ; j++ ) {
		// 		if (  j%50 == 0 && j != 0 )
		// 		  puts("");
		// 		printf("%c",seq[j]);
  //     }
  //     puts("");

  //     printf(">%sp2\n",name);
  //     for ( int j = 0 ; j < reallen ; j++ ) {
		// 		if (  j%50 == 0 && j != 0 )
	 //  			puts("");
		// 		printf("%c",seq[j]);
  //     }
  //     puts("");
  //   }
    
  //   free ( seq );
  // }

  
  
  // fprintf(save,"chr\tchrStart\tchrEnd\tidOfDupl\tnewChr\tindexInNewChr\tstrand\tisReversed\n");
  
  // vector < pss > dupl = read_tab ();
  // map < int , bool > used;
  // vector < char > rnd1, rnd2;

  // for ( int i = 0 ; i < 10000 ; i++ ) {
  //   rnd1.pb ( 'N' );
  //   rnd2.pb ( 'N' );;
  // }

  int len_of_first, len_of_second;
  char *seq = faidx_fetch_seq ( ref_fai, seg1.chr.c_str(), seg1.beg + aft, seg1.beg + aft+len-1, &len_of_first );
  char *seq2 = faidx_fetch_seq ( ref_fai, seg2.chr.c_str(), seg2.beg, seg2.beg+len-1, &len_of_second );

  if ( len_of_first != len_of_second )
  	uneven_cnt++;


  // STRAAAAAAAAAAAANNNNNNNNNDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
  for ( int i = 0 ; i < min (len_of_first, len_of_second) ; i++ )
 		if ( seq[i] == seq2[i] ) { // STRANDLA ALAKALI BISILER EKSIKKKKKKKK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 			set < segment > :: iterator it;
 			segment tmp ( seg1.chr, seg1.beg + aft + i , seg1.beg + aft + i, 0 );
 			it = suns.find ( tmp );
 			if ( it != suns.end() )
 				suns.erase (it);

 			segment tmp2 ( seg2.chr, seg2.beg + i , seg2.beg + i, 0 );
 			it = suns.find ( tmp2 );
 			if ( it != suns.end() )
 				suns.erase (it);
 		}

  free (seq);
  free (seq2);

  // const int MAXX = (int) 2e9;
  /*
  for ( int i = 0 ; i < dupl.size() ; i++ ) {

    if ( rnd1.size() > MAXX || rnd2.size() > MAXX ) break;
    if ( dupl[i][0].size() > 5 || dupl[i][6].size() > 5 ) continue;
    
    vector < string > cur = dupl[i];
    int len, len2, len3, id;
    
    id = atoi ( cur[10].c_str() );
    if ( used[id] == true ) continue;
    used[id] = true;
    
    char *seq = faidx_fetch_seq ( ref_fai , cur[0].c_str() , atoi(cur[1].c_str()), atoi(cur[2].c_str()),  &len );
    char *seq2 = faidx_fetch_seq ( ref_fai , cur[6].c_str() , atoi(cur[7].c_str()), atoi(cur[8].c_str()),  &len2 );
    char *seq3 = faidx_fetch_seq ( ref_fai , cur[0].c_str() , atoi(cur[1].c_str()), atoi(cur[2].c_str()),  &len3 );

    for ( int j = 0 ; j < len3 ; j++ )
      seq3[j] = toupper ( seq3[j] );
    for ( int j = 0 ; j < len3/2 ; j++ ) {
      char tmpp = seq3[j];
      seq3[j] = seq3[len3-j-1];
      seq3[len3-j-1] = tmpp;
    }
    for ( int j = 0 ; j < len3 ; j++ ) {
      if ( seq3[j] == 'N' ) continue;
      if ( seq3[j] == 'A' ) seq3[j] = 'T';
      else if ( seq3[j] == 'T' ) seq3[j] = 'A';
      else if ( seq3[j] == 'G' ) seq3[j] = 'C';
      else if ( seq3[j] == 'C' ) seq3[j] = 'G';
    }

    int chc = rand()%6;
    if ( chc == 5 ) {
      if ( rnd1.size() <= rnd2.size() ) {
	fprintf(save,"%s\t%d\t%d\t%d\tchrRn1\t%d\t%s\t+\n",cur[0].c_str(),atoi(cur[1].c_str()), atoi(cur[2].c_str()),id,(int)rnd1.size(),cur[5].c_str());
	for ( int j = 0 ; j < len3 ; j++ )
	  rnd1.pb ( seq3[j] );
      }
      else {
	fprintf(save,"%s\t%d\t%d\t%d\tchrRn2\t%d\t%s\t+\n",cur[0].c_str(),atoi(cur[1].c_str()), atoi(cur[2].c_str()),id,(int)rnd2.size(),cur[5].c_str());
	for ( int j = 0 ; j < len3 ; j++ )
	  rnd2.pb ( seq3[j] );
      }
    }
    if ( chc == 0  || chc == 3 || chc == 4 ) {
      fprintf(save,"%s\t%d\t%d\t%d\tchrRn1\t%d\t%s\t-\n",cur[0].c_str(),atoi(cur[1].c_str()), atoi(cur[2].c_str()),id,(int)rnd1.size(),cur[5].c_str());
      fprintf(save,"%s\t%d\t%d\t%d\tchrRn2\t%d\t%s\t-\n",cur[6].c_str(),atoi(cur[7].c_str()), atoi(cur[8].c_str()),id,(int)rnd2.size(),cur[5].c_str());
      for ( int j = 0 ; j < len ; j++ )
	rnd1.pb ( seq[j] );
      for ( int j = 0 ; j < len2 ; j++ )
	rnd2.pb ( seq2[j] );
    }
    if ( chc == 1 ) {
      if ( rnd1.size() <= rnd2.size() ) {
	fprintf(save,"%s\t%d\t%d\t%d\tchrRn1\t%d\t%s\t-\n",cur[0].c_str(),atoi(cur[1].c_str()), atoi(cur[2].c_str()),id,(int)rnd1.size(),cur[5].c_str());
	for ( int j = 0 ; j < len ; j++ )
	  rnd1.pb ( seq[j] );
      }
      else {
	fprintf(save,"%s\t%d\t%d\t%d\tchrRn2\t%d\t%s\t-\n",cur[0].c_str(),atoi(cur[1].c_str()), atoi(cur[2].c_str()),id,(int)rnd2.size(),cur[5].c_str());
	for ( int j = 0 ; j < len ; j++ )
	  rnd2.pb ( seq[j] );
      }
    }
    if ( chc == 2 ) {
      if ( rnd2.size() <= rnd1.size() ) {
	fprintf(save,"%s\t%d\t%d\t%d\tchrRn2\t%d\t%s\t-\n",cur[6].c_str(),atoi(cur[7].c_str()), atoi(cur[8].c_str()),id,(int)rnd2.size(),cur[5].c_str());
	for ( int j = 0 ; j < len2 ; j++ )
	  rnd2.pb ( seq2[j] );
      }
      else {
	fprintf(save,"%s\t%d\t%d\t%d\tchrRn1\t%d\t%s\t-\n",cur[6].c_str(),atoi(cur[7].c_str()), atoi(cur[8].c_str()),id,(int)rnd1.size(),cur[5].c_str());
	for ( int j = 0 ; j < len2 ; j++ )
	  rnd1.pb ( seq2[j] );
      }
    }
    
    free ( seq );
    free ( seq2 );
    free ( seq3 );
  }
  */
  
  
}

segment current_grp ( segment cur ) {
	if ( grp[cur] == cur )
		return cur;
	return grp[cur] = current_grp ( grp[cur] );
}

void find_groups () {
	for ( int i = 0 ; i < pairs.size() ; i++ ) {
		if ( current_grp ( pairs[i].first ) == current_grp ( pairs[i].second ) )
			continue;
		grp[ grp[pairs[i].second] ] = grp[pairs[i].first];
	}

	for ( auto i: grp ) {
		super_segments.insert ( i.second );
		group_members[i.second].insert ( i.first );
	}
}

void find_variations () {
	set < segment > :: iterator it, it2;
	for ( it = all_segments.begin() ; it != all_segments.end() ; it++ ) {
		it2 = it;
		it2++;
		for ( ; it2 != all_segments.end() ; it++ ) {
			if ( it2->chr.compare(it->chr) != 0 || (it2->beg) > (it->en) )
				break;
			int aft = (it2->beg) - (it->beg);
			int lenn = min( it2->en, it->en ) - (it2->beg) + 1;
			similar tmp;
			tmp.seg1 = *it;
			tmp.seg2 = *it2;
			tmp.after_first = aft;
			tmp.len = lenn;
			similars.push_back ( tmp );
		}
	}
}

void find_suns () {

	fprintf(stderr,"Before deleting any, sun size = %d\n",(int)suns.size());

	for ( auto i: group_members ) {
		set < segment > &cur_grp = i.second;
		set < segment > :: iterator it, it2;

		for ( it = cur_grp.begin() ; it != cur_grp.end() ; it++ ) {
			it2 = it;
			it2++;
			while ( it2 != cur_grp.end() ) {
				compare ( *it, *it2, 0, min ( (it->en) - (it->beg) + 1, (it2->en) - (it2->beg) + 1 ) );
				it2++;
			}
		}
	}

	fprintf(stderr,"After deleting inter-groups, sun size = %d\n",(int)suns.size());

	for ( int i = 0 ; i < similars.size () ; i++ ) {
		similar &tmp = similars[i];

		segment sup1 = current_grp (tmp.seg1), sup2 = current_grp (tmp.seg2);
		set < segment > &grp1 = group_members[sup1];
		set < segment > &grp2 = group_members[sup2];
		set < segment > :: iterator it, it2;

		for ( it = grp1.begin() ; it != grp1.end() ; it++ ) {
			for ( it2 = grp2.begin() ; it2 != grp2.end() ; it2++ )
				compare ( *it, *it2, tmp.after_first, tmp.len );
		}

	}

	fprintf(stderr,"After deleting variations, sun size = %d\n",(int)suns.size());
	fprintf(stderr,"fai_fetch uneven length count = %d\n",uneven_cnt);
}

int main( int argc, char** argv ) {

  srand( time (NULL) );
  get_file_name ( argc, argv );
  pairs = read_tab();

  find_groups();
  find_variations();

  ref_fai = fai_load (files.fasta_file);
  save = fopen(files.output_file,"w");
  find_suns ();
  fai_destroy ( ref_fai );

  // STRAND HALLEDIP, PRINTLEME KALDI

  fclose( save );
  return 0;
}

vector < pss > read_tab () {

  vector < pss > ret;
  ifstream file(files.duplication_file);
  
  if (file.is_open()) {
    
    string line;
    while (getline(file, line)) {

      vector < string > words;
      int st,en,i = -1;

      while ( i != line.size() ) {
				string tmp;
				for ( i++ ; i != line.size() ; i++ ) {
				  if ( line[i] == '\t' || line[i] == ' ' || line[i] == '\n' ) break;
				  tmp += line[i];
				}
				words.pb ( tmp );
      }

      int strand1 = 0, strand2 = 0;
      if ( words[8].compare("+") == 0 )
				strand1 = 1;
      if ( words[8].compare("-") == 0 )
				strand1 = 2;

      if ( words[9].compare("+") == 0 )
				strand2 = 1;
      if ( words[9].compare("-") == 0 )
				strand2 = 2;
      segment tmp1 (words[0],stoi(words[1]),stoi(words[2]),strand1);
      segment tmp2 (words[3],stoi(words[4]),stoi(words[5]),strand2);
      ret.pb ( make_pair(tmp1, tmp2) );
      all_segments.insert (tmp1);
      all_segments.insert (tmp2);
      grp[tmp1] = tmp1;
      grp[tmp2] = tmp2;

      for ( int j = stoi(words[1]) ; j <= stoi(words[2]) ; j++ ) {
      	segment nucl ( words[0], j, j, 0 );
      	suns.insert ( nucl );
      }
      for ( int j = stoi(words[4]) ; j <= stoi(words[5]) ; j++ ) {
      	segment nucl ( words[3], j, j, 0 );
      	suns.insert ( nucl );
      }
    }
    
    file.close();
  }
  else {
    fprintf(stderr,"Unable to open file %s. Make sure the system has read permition on this file and run the program again.\n",files.duplication_file);
    exit(0);
  }

  return ret;
}
