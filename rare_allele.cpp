#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <list>
#include <algorithm>


using namespace std;


typedef unsigned char uchar;
using vvchar = vector<vector<char> >;

vector<string> pullNames(string);
uchar ** convertVecToData(vvchar &);
double getGenCol(string line, int col, double & phys);

bool nearZero(double a) {
	if(a <= 0.000001 && a >= -0.000001) { return 1; }
	return 0;
}

vector<string> split(string s, string delim) {
	vector<string> out;
	int pos = 0;
	string token;
	while ((pos = s.find(delim)) != string::npos) {
		token = s.substr(0,pos);
		out.push_back(token);
		s.erase(0,pos+delim.length());
	}
	out.push_back(s);
	return out;
}

void addToVector(vector<string> & a, string b) {
	for(int i = 0; i < a.size(); i++) {
		if(a[i].compare(b) == 0) { return; }
	}
	a.push_back(b);

}

string getCategory(string s) {
	if(s.find("HIGH") != string::npos) {
		return "HIGH";
	} else if(s.find("MODERATE") != string::npos) {
		return "MODERATE";
	} else if(s.find("LOW") != string::npos) {
		return "LOW";
	}
	return "MODIFIER";
}

string typeLabel(vector<string> & info_s) {
	string out;
	vector<string> vout;
	string category;
	for(int i = 0; i < info_s.size(); i++) {
		vector<string> s = split(info_s[i],"=");
		if(s[0].compare("ANN") == 0) {
			category = getCategory(s[1]);
			vector<string> anns = split(s[1],",");
			for(int j = 0; j < anns.size(); j++) {
				vector<string> al = split(anns[j],"|");
				addToVector(vout,al[1]);
			}
		}
	}
	out += (category + ",");
	for(int i = 0; i < vout.size(); i++) {
		out += vout[i];
		if(i != vout.size() - 1) { out += ","; }
	}
	return out;


}

bool isCpg(vector<string> & info_s) {
	for(int i = 0; i < info_s.size(); i++) {
		if(info_s[i].compare("CPG_TAG") == 0) { return 1; }
	}
	return 0;

}

bool isIndel(vector<string> & info_s) {
	for(int i = 0; i < info_s.size(); i++) {
		if(info_s[i].compare("INDEL") == 0) { return 1; }
	}
	return 0;
}



class maxhapdata {
	public:
	vector<int> cutoffs;
	int maxhap;
	maxhapdata() {

	};
	maxhapdata(int size) {
		vector<int> t(size,-1);
		cutoffs = t;
		maxhap = -1;
	}
	
};

class hapset {
	public:
	int hap_id; //Even number, hap_id and hap_id+1 are idxs of indiv with rare allele
	int snp_id;
	list<int> hap_idx_list;
	maxhapdata data_left;
	maxhapdata data_right;
	hapset(int hap, int snp, int hap_count, list<int> & hil) {
		hap_id = hap;
		snp_id = snp;
		hap_idx_list = hil;
		maxhapdata t1(hap_count);
		maxhapdata t2(hap_count);
		data_left = t1;
		data_right = t2;
	}
	hapset(int hap, int snp, int hap_count, list<int> & hil, int extra) {
		hap_id = hap;
		snp_id = snp;
		hap_idx_list = hil;
		hap_idx_list.push_back(extra);
		hap_idx_list.sort();
		maxhapdata t1(hap_count);
		maxhapdata t2(hap_count);
		data_left = t1;
		data_right = t2;
	}
	void fillLeft(int mn) {
		for(auto i = hap_idx_list.begin(); i != hap_idx_list.end(); i++) {
			if(data_left.cutoffs[*i] == -1) {
				data_left.cutoffs[*i] = mn;
			}
		}
	}
	void fillRight(int mx) {
		for(auto i = hap_idx_list.begin(); i != hap_idx_list.end(); i++) {
			if(data_right.cutoffs[*i] == -1) {
				data_right.cutoffs[*i] = mx;
			}
		}
	}
	void print() {
		cout << "Hap id: " << hap_id << endl;
		cout << "Snp id: " << snp_id << endl;
		cout << "Max: " << data_left.maxhap << " " << data_right.maxhap << endl;
		list<int>::iterator itr = hap_idx_list.begin();
		for(int i = 0; i < 10; i++) {
			cout << *itr << " " << data_left.cutoffs[i] << " " << data_right.cutoffs[i] << endl;
			++itr;
		}
	}

};

class snp_data {
	public: 
	int allele_count;
	int position;
	int phased;
	double genetic_position;
	string annotation;
	snp_data(int ac, int pos, int phs) {
		allele_count = ac;
		position = pos;
		phased = phs;
	}
	snp_data(int ac, int pos, int phs, string ann) {
		allele_count = ac;
		position = pos;
		phased = phs;
		annotation = ann;
	}
	bool isSingleton() {
		if(allele_count == 1) {
			return 1;
		}
		return 0;
	}
	bool isValidForComp() {
		if(allele_count > 2) {
			return 1;
		}
		//Add check for CpG flag
		return 0;
	}
};

class vcf_data {
	public:
	uchar ** data;
	int hap_count;
	int snp_count;
	int transpose;
	int compress;
	int gendist_minmax[2];
	vector<snp_data> snps;
	vector<string> indiv_list;
	char getHap(int hap, int snp) {
		return char(data[hap][snp]);
	}
	void readVcf(string filename) {
		ifstream infile;
		infile.open(filename.c_str());
		string line;
		vvchar hap_hold;
		bool vector_set = 0;
		while(getline(infile,line)) {
			if(line[0] == '#') {
				if(line[1] != '#') {
					indiv_list = pullNames(line);
				}
				continue;
			}
			if(line.size() == 0) {
				break;
			}
			stringstream s(line);
			string junk, ref, alt, hap, info;
			int phys;
			s >> junk >> phys >> junk >> ref >> alt >> junk >> junk >> info;
			vector<string> info_s = split(info,";");
			if(ref.size() != 1 || isIndel(info_s) || isCpg(info_s)) {
				continue;
			}
			int ac = 0;
			int phased = -1;
			int ind_idx = 0;
			string annotation = typeLabel(info_s);
			s >> junk;
			while(s >> hap) {
				if(vector_set == 0) {
					vector<char> t;
					hap_hold.push_back(t);
					vector<char> t2;
					hap_hold.push_back(t2);
				}
				if(phased == -1) { phased = (hap[1] == '|') ? 1 : 0; }
				if(hap[0] == '0') { 
					hap_hold[ind_idx++].push_back(0);
				} else {
					hap_hold[ind_idx++].push_back(1);
					ac++;
				}
				if(hap[2] == '0') {
					hap_hold[ind_idx++].push_back(0);
				} else {
					hap_hold[ind_idx++].push_back(1);
					ac++;
				}
			}
			snp_data snp(ac,phys,phased,annotation);
			snps.push_back(snp);
			vector_set = 1;
		}
		snp_count = snps.size();
		hap_count = hap_hold.size();
		data = convertVecToData(hap_hold);
	}
	void printVcf() {
		cout << "Hap count: " << hap_count << endl;
		cout << "Snp count: " << snp_count << endl;
		for(int i = 0; i < snp_count; i++) {
			cout << "Snp " << i << ", position " << snps[i].position << ", genetic " << snps[i].genetic_position << " " << gendist_minmax[0] << " " << gendist_minmax[1] << " " ;
			for(int j = 0; j < hap_count; j++) {
				cout << int(data[j][i]);
			}
			cout << endl;
		}
	}
	int getSnpDist(int base, int end) {
		return snps[base].position - snps[end].position;
	}
	double getSnpGm(int base, int end) {
		return snps[base].genetic_position - snps[end].genetic_position;
	}
	void addGenData(string filename, int offset) {
		ifstream gm;
		gm.open(filename.c_str());
		string line;
		vector<double> g0,g1;
		double a,b;
		getline(gm,line);
		while(gm) {
			getline(gm,line);
			b = getGenCol(line,offset,a);
			if(gm.eof()) { break; }
			if(b <= -0.5) {
				cerr << "Error with genetic map\n";
				cerr << a << " " << b << " " << line << endl;
				exit(0);
			}
			//if(g1.size() == 0 || b != g1[g1.size()-1]) {
			if(!nearZero(b) && (g1.size() < 2 || b != g1[g1.size()-1])) {
				g0.push_back(a);
				g1.push_back(b);
				//cout << a << " " << b << endl;
			}
		}
		cerr << "Finding min\n";
		for(int ii = 1; ii < g0.size(); ii++) {
			if(g1[ii] != 0 && g1[ii-1] == 0) {
				gendist_minmax[0] = g0[ii];
				break;
			}
		}
		cerr << "Finding max\n";
		for(int ii = g0.size() - 1; ii > 0; ii--) {
			if(g1[ii] != g1[ii-1]) {
				gendist_minmax[1] = g0[ii];
			}
		}
		cerr << "Making distances\n";
		for(int i = 0; i < snp_count; i++) {
			int j;
			if(snps[i].position <= g0[0]) {
				double startrate = (g1[1]-g1[0])/(g0[1]-g0[0]);
				snps[i].genetic_position = startrate*(snps[i].position-g0[0]);
			} else {
				for(j = 0; j < g0.size()-1; j++) {
					if(snps[i].position <= g0[j]) {
						snps[i].genetic_position = g1[0];
						break;
					}
					if((snps[i].position >= g0[j]) && (snps[i].position <= g0[j+1])) {
						snps[i].genetic_position = g1[j]+(g1[j+1]-g1[j])/(g0[j+1]-g0[j])*(snps[i].position-g0[j]);
						break;
					}
				}
				if(j == g0.size() - 1) {
					double endrate = (g1[g1.size()-1]-g1[g1.size()-2])/(g0[g0.size()-1]-g0[g0.size()-2]);
					snps[i].genetic_position = snps[i-1].genetic_position+endrate*(snps[i].position - snps[i-1].position);
				}
			}
		}
	}
};

string ast(int a, bool b) {
	//cout << a;
	stringstream s("");
	s << a;
	if(b) {
		s << '*';
	}
	s << '\t';
	return s.str();
}

string ast(double a, bool b) {
	//cout << a;
	stringstream s("");
	s << a;
	if(b) {
		s << '*';
	}
	s << '\t';
	return s.str();
}


class output_data {
	public:
	int pos;
	int p1_5_bp, p2_5_bp, p1_3_bp, p2_3_bp, p1_total_bp, p2_total_bp;
	double p1_5_gen, p2_5_gen, p1_3_gen, p2_3_gen, p1_total_gen, p2_total_gen;
	bool same_1,same_2;
	bool end_flags[12] = {false};
	int single_phase;
	string ann;
	output_data(vcf_data & vcf, hapset & h1, hapset & h2) {
		int front_end = 0;
		int back_end = vcf.snp_count-1;
		single_phase = 0;
		ann = vcf.snps[h1.snp_id].annotation;
		pos = vcf.snps[h1.snp_id].position;
		p1_5_bp = vcf.getSnpDist(h1.data_right.maxhap,h1.snp_id);
		if(h1.data_right.maxhap == back_end) { end_flags[0] = 1; end_flags[4] = 1; }
		p1_3_bp = vcf.getSnpDist(h1.snp_id,h1.data_left.maxhap);
		if(h1.data_left.maxhap == front_end) { end_flags[1] = 1; end_flags[5] = 1; }
		p2_5_bp = vcf.getSnpDist(h2.data_right.maxhap,h2.snp_id);
		if(h2.data_right.maxhap == back_end) { end_flags[2] = 1; end_flags[6] = 1; }
		p2_3_bp = vcf.getSnpDist(h2.snp_id,h2.data_left.maxhap);
		if(h2.data_left.maxhap == front_end) { end_flags[3] = 1; end_flags[7] = 1; }
		p1_5_gen = vcf.getSnpGm(h1.data_right.maxhap,h1.snp_id);
		p1_3_gen = vcf.getSnpGm(h1.snp_id,h1.data_left.maxhap);
		p2_5_gen = vcf.getSnpGm(h2.data_right.maxhap,h2.snp_id);
		p2_3_gen = vcf.getSnpGm(h2.snp_id,h2.data_left.maxhap);
		p1_total_bp = 0, p2_total_bp = 0, p1_total_gen = 0, p2_total_gen = 0;
		int phys_1_idx = 0, phys_2_idx = 0, gen_1_idx = 0, gen_2_idx = 0;
		for(auto itr = h1.hap_idx_list.begin(); itr != h1.hap_idx_list.end(); ++itr) {
			int i = *itr;
			int t = vcf.getSnpDist(h1.data_right.cutoffs[i],h1.data_left.cutoffs[i]);
			if (t > p1_total_bp) {
				p1_total_bp = t;
				phys_1_idx = i;
			}
			/*t = vcf.getSnpDist(h2.data_right.cutoffs[i],h2.data_left.cutoffs[i]);
			if (t > p2_total_bp) {
				p2_total_bp = t;
				phys_2_idx = i;
			}*/
			double g = vcf.getSnpGm(h1.data_right.cutoffs[i],h1.data_left.cutoffs[i]);
			if (g > p1_total_gen) {
				p1_total_gen = g;
				gen_1_idx = i;
			}
			/*g = vcf.getSnpGm(h2.data_right.cutoffs[i],h2.data_left.cutoffs[i]);
			if (g > p2_total_gen) {
				p2_total_gen = g;
				gen_2_idx = i;
			}*/
			if(h1.data_right.cutoffs[phys_1_idx] == back_end || h1.data_left.cutoffs[phys_1_idx] == front_end) { end_flags[8] = 1; }
			//if(h2.data_right.cutoffs[phys_2_idx] == back_end || h2.data_left.cutoffs[phys_2_idx] == front_end) { end_flags[9] = 1; }
			if(h1.data_right.cutoffs[gen_1_idx] == back_end || h1.data_left.cutoffs[gen_1_idx] == front_end) { end_flags[10] = 1; }
			//if(h2.data_right.cutoffs[gen_2_idx] == back_end || h2.data_left.cutoffs[gen_2_idx] == front_end) { end_flags[11] = 1; }
		}
		for(auto itr = h2.hap_idx_list.begin(); itr != h2.hap_idx_list.end(); ++itr) {
			int i = *itr;
			int t = vcf.getSnpDist(h2.data_right.cutoffs[i],h2.data_left.cutoffs[i]);
			if (t > p2_total_bp) {
				p2_total_bp = t;
				phys_2_idx = i;
			}
			double g = vcf.getSnpGm(h2.data_right.cutoffs[i],h2.data_left.cutoffs[i]);
			if (g > p2_total_gen) {
				p2_total_gen = g;
				gen_2_idx = i;
			}
			if(h2.data_right.cutoffs[phys_2_idx] == back_end || h2.data_left.cutoffs[phys_2_idx] == front_end) { end_flags[9] = 1; }
			if(h2.data_right.cutoffs[gen_2_idx] == back_end || h2.data_left.cutoffs[gen_2_idx] == front_end) { end_flags[11] = 1; }
		}
		same_1 = (phys_1_idx == gen_1_idx);
		same_2 = (phys_2_idx == gen_2_idx);

	}
	output_data(vcf_data & vcf, hapset & h1) {
		int front_end = 0;
		int back_end = vcf.snp_count-1;
		single_phase = 1;
		ann = vcf.snps[h1.snp_id].annotation;
		pos = vcf.snps[h1.snp_id].position;
		p1_5_bp = vcf.getSnpDist(h1.data_right.maxhap,h1.snp_id);
		if(h1.data_right.maxhap == back_end) { end_flags[0] = 1; end_flags[4] = 1; }
		p1_3_bp = vcf.getSnpDist(h1.snp_id,h1.data_left.maxhap);
		if(h1.data_left.maxhap == front_end) { end_flags[1] = 1; end_flags[5] = 1; }
		p2_5_bp = 0;
		p2_3_bp = 0;
		p1_5_gen = vcf.getSnpGm(h1.data_right.maxhap,h1.snp_id);
		p1_3_gen = vcf.getSnpGm(h1.snp_id,h1.data_left.maxhap);
		p2_5_gen = 0;
		p2_3_gen = 0;
		p1_total_bp = 0, p2_total_bp = 0, p1_total_gen = 0, p2_total_gen = 0;
		int phys_1_idx = 0, phys_2_idx = 0, gen_1_idx = 0, gen_2_idx = 0;
		for(auto itr = h1.hap_idx_list.begin(); itr != h1.hap_idx_list.end(); ++itr) {
			int i = *itr;
			int t = vcf.getSnpDist(h1.data_right.cutoffs[i],h1.data_left.cutoffs[i]);
			if (t > p1_total_bp) {
				p1_total_bp = t;
				phys_1_idx = i;
			}
			p2_total_bp = 0;
			double g = vcf.getSnpGm(h1.data_right.cutoffs[i],h1.data_left.cutoffs[i]);
			if (g > p1_total_gen) {
				p1_total_gen = g;
				gen_1_idx = i;
			}
			p2_total_gen = 0;
			if(h1.data_right.cutoffs[phys_1_idx] == back_end || h1.data_left.cutoffs[phys_1_idx] == front_end) { end_flags[8] = 1; }
			if(h1.data_right.cutoffs[gen_1_idx] == back_end || h1.data_left.cutoffs[gen_1_idx] == front_end) { end_flags[10] = 1; }
		}
		same_1 = (phys_1_idx == gen_1_idx);
		same_2 = 1;

	}
	output_data() { }
	void printData() {
		//cout << pos << '\t' << p1_5_bp << '\t' << p2_5_bp << '\t' << p1_3_bp << '\t' << p2_3_bp << '\t' << p1_5_gen << '\t' << p2_5_gen << '\t' << p1_3_gen << '\t' << p2_3_gen << '\t' << p1_total_bp << '\t' << p2_total_bp << '\t' << p1_total_gen << '\t' << p2_total_gen << '\t' << same_1 << '\t' << same_2 << endl;
		cout << pos << '\t' << ast(p1_5_bp,end_flags[0]) << ast(p1_3_bp,end_flags[1]) << ast(p2_5_bp,end_flags[2]) << ast(p2_3_bp,end_flags[3]) << ast(p1_5_gen,end_flags[4]) << ast(p1_3_gen,end_flags[5]) << ast(p2_5_gen,end_flags[6]) << ast(p2_3_gen,end_flags[7]) << ast(p1_total_bp,end_flags[8]) << ast(p2_total_bp,end_flags[9]) << ast(p1_total_gen,end_flags[10]) << ast(p2_total_gen,end_flags[11]) << same_1 << '\t' << same_2 << '\t' << ann << endl;
	}



};

uchar ** convertVecToData(vvchar & hap_hold) {
	int hap_count = hap_hold.size();
	int snp_count = hap_hold[0].size();
	uchar ** d = new uchar*[hap_count];
	for(int i = 0; i < hap_hold.size(); i++) {
		d[i] = new uchar[snp_count];
		for(int j = 0; j < snp_count; j++) {
			d[i][j] = hap_hold[i][j];
		}
		hap_hold[i].clear();
		(vector<char>()).swap(hap_hold[i]);
	}
	return d;


}




vector<string> pullNames(string line) {
        stringstream s(line);
        string hold, junk;
        s >> junk >> junk >> junk >> junk >> junk >> junk >> junk >> junk >> junk; //9 junks, verify that merged files didn't skip a person
        vector<string> names;
        while(s >> hold) {
                names.push_back(hold);
        }
        return names;
}

double getGenCol(string line, int col, double & phys) {
	stringstream s(line);
	s >> phys;
	double j = 0, out = 0;
	for(int i = 0; i < col; i++) {
		s >> j;
	}
	bool b = (s >> out);
	if(!b) { return -1.0; }
	return out;
}

list<int> subsampleIdx(vector<int> & vals, int subsample) {
	random_shuffle(vals.begin(), vals.end());
	list<int> out;
	for(int i = 0; i < subsample; i++) {
		out.push_back(vals[i]);
	}
	out.sort();
	return out;
}


list<int> getSampleListSingle(vcf_data & vcf, int idx, vector<int> & rare_idx, int subsample) {
	list<int> out; 
	vector<int> vhold;
	for(int i = 0; i < vcf.hap_count; i+=2) {
		char h1 = vcf.getHap(i,idx);
		char h2 = vcf.getHap(i+1,idx);
		if(h1 == 0 && h2 == 0) {
			vhold.push_back(i);
			vhold.push_back(i+1);
		} else {
			rare_idx.push_back(i);
			rare_idx.push_back(i+1);
		}	
	}
	if(subsample == 0) {
		//cout << "Adding to list\n";
		for(int i = 0; i < vhold.size(); i++) {
			out.push_back(vhold[i]);
		}
	} else {
		out = subsampleIdx(vhold,subsample);
	}
	return out;
}

//For subsampling, 0 = dont, anything positive means list will contain subsample number of indices
list<int> getSampleListMult(vcf_data & vcf, int idx, vector<int> & rare_single, vector<int> & rare_double, int subsample) {
	list<int> out;
	vector<int> vhold;
	for(int i = 0; i < vcf.hap_count; i+=2) {
		char h1 = vcf.getHap(i,idx);
		char h2 = vcf.getHap(i+1,idx);
		if(h1 == 0 && h2 == 0) {
			vhold.push_back(i);
			vhold.push_back(i+1);
		} else if(h1 == 1 && h2 == 1) {
			rare_double.push_back(i);
			rare_double.push_back(i+1);
		} else if(h1 == 1) {
			rare_single.push_back(i);
			vhold.push_back(i+1);
		} else {
			rare_single.push_back(i+1);
			vhold.push_back(i);
		}
		
	}
	//Subsample here
	if(subsample == 0) {
		//cout << "Adding to list\n";
		for(int i = 0; i < vhold.size(); i++) {
			out.push_back(vhold[i]);
		}
	} else {
		out = subsampleIdx(vhold,subsample);
	}

	return out;
}

list<int> getModList(list<int> & haps, int sw) {
	list<int> out;
	for(auto i = haps.begin(); i != haps.end(); ++i) {
		if(*i != sw) { out.push_back(*i); }
		else {
			if(sw%2 == 1) {
				out.push_back(*i-1);
			} else {
				out.push_back(*i+1);
			}
			
		}
	}
	return out;

}

bool fillStoppingPoint(vcf_data & vcf, hapset & h, maxhapdata & mhd, list<int> & comp_haps, int snp_idx) {
	list<int>::iterator i = comp_haps.begin();
	while(i != comp_haps.end()) {
		char rare_hap = vcf.getHap(h.hap_id,snp_idx);
		char current_hap = vcf.getHap(*i,snp_idx);
		if(rare_hap != current_hap) {
			mhd.cutoffs[*i] = snp_idx;
			comp_haps.erase(i++);
		} else { 
			++i;
		}

	}
	if(comp_haps.size() == 0) { return 1; } //Return 1 means loop in calling function should end
	return 0;
}


void findMax(vcf_data & vcf, hapset & h) {
	list<int> c_left = h.hap_idx_list;
	//left portion
	int i = h.snp_id-1;
	for(i = h.snp_id-1; i > 0; i--) {
		if(!vcf.snps[i].isValidForComp()) { continue; }
		bool filled = fillStoppingPoint(vcf,h,h.data_left,c_left,i);
		if(filled) {
			break;
		}
	}
	h.data_left.maxhap = i;
	//cout << i << endl;
	if(i == 0) {
		h.fillLeft(0);
		//cout << "filling left side\n";
	}
	list<int> c_right = h.hap_idx_list;
	i = h.snp_id+1;
	for(i = h.snp_id+1; i < vcf.snp_count - 1; i++) {
		if(!vcf.snps[i].isValidForComp()) { continue; }
		bool filled = fillStoppingPoint(vcf,h,h.data_right,c_right,i);
		if(filled) {
			break;
		}
	}
	h.data_right.maxhap = i;
	//cout << i << " " << vcf.snp_count << endl;
	if(i >= vcf.snp_count-1) {
		h.fillRight(vcf.snp_count-1);
		h.data_right.maxhap = vcf.snp_count-1;
		//cout << "filling right side " << vcf.snp_count << endl;
	}

}

void printVerify(vcf_data & vcf, hapset & h) {
	cout << h.snp_id << endl;
	cout << h.hap_id << '\t';
	for(int i = h.snp_id; i >= h.data_left.maxhap; i--) {
		cout << int(vcf.getHap(h.hap_id,i));
	}
	cout << endl;
	for(auto itr = h.hap_idx_list.begin(); itr != h.hap_idx_list.end(); ++itr) {
		cout << *itr << '\t';
		for(int i = h.snp_id; i >= h.data_left.cutoffs[*itr]; i--) {
			cout << int(vcf.getHap(*itr,i));
		}
		cout << endl;
	}
	cout << endl;
	cout << h.hap_id << '\t';
	for(int i = h.snp_id; i <= h.data_right.maxhap; i++) {
		cout << int(vcf.getHap(h.hap_id,i));
	}
	cout << endl;
	for(auto itr = h.hap_idx_list.begin(); itr != h.hap_idx_list.end(); ++itr) {
		cout << *itr << '\t';
		for(int i = h.snp_id; i <= h.data_right.cutoffs[*itr]; i++) {
			cout << int(vcf.getHap(*itr,i));
		}
		cout << endl;
	}
}


void findLongestHaps(vcf_data & vcf, int idx) {
	vector<int> rare_idx;
	vector<hapset> vh;
	if(vcf.snps[idx].allele_count == 1) {
		//cout << "SINGLETON\n";
		vector<int> rare_idx;
		list<int> comparison_haps = getSampleListSingle(vcf,idx,rare_idx,0);
		//cout << "SINGLE " << rare_idx[0] << " " << rare_idx[1] << endl;
		hapset h1(rare_idx[0],idx,vcf.hap_count,comparison_haps,rare_idx[1]);
		hapset h2(rare_idx[1],idx,vcf.hap_count,comparison_haps,rare_idx[0]);
		findMax(vcf,h1);
		findMax(vcf,h2);
		output_data od(vcf,h1,h2);
		od.printData();
	} else { 
		//cout << "MULTI-COPY " << vcf.snps[idx].allele_count << endl;
		vector<int> rare_single, rare_double;
		list<int> comparison_haps = getSampleListMult(vcf,idx,rare_single,rare_double,0);
		//for(int ip = 0; ip < rare_single.size(); ip++) {
		//	cout << rare_single[ip] << endl;
		//}
		for(int i = 0; i < rare_single.size(); i++) {
			hapset h1(rare_single[i],idx,vcf.hap_count,comparison_haps);
			int sw_person;
			if(rare_single[i]%2==1) {
				sw_person = rare_single[i]-1;
			} else {
				sw_person = rare_single[i]+1;
			}
			list<int> swap_phase_haps = getModList(comparison_haps,sw_person);
			hapset h2(sw_person,idx,vcf.hap_count,swap_phase_haps);
			findMax(vcf,h1);
			findMax(vcf,h2);
			output_data od(vcf,h1,h2);
			od.printData();
		}
		for(int i = 0; i < rare_double.size(); i+=2) {
			hapset h1(rare_double[i],idx,vcf.hap_count,comparison_haps);
			findMax(vcf,h1);
			output_data od(vcf,h1);
			od.printData();
			hapset h2(rare_double[i+1],idx,vcf.hap_count,comparison_haps);
			findMax(vcf,h2);
			output_data od2(vcf,h2);
			od2.printData();
		}
	}

}
		




int main(int argc, char ** argv) {
	setbuf(stdout,NULL);
	srand(1);
	string vcf_name = argv[1];
	vcf_data v;
	v.readVcf(vcf_name);
	string map_name = argv[2];
	v.addGenData(map_name,0);
	//v.printVcf();
	for(int i = 0; i < v.snps.size(); i++) {
		if(v.snps[i].allele_count <= 5 && v.snps[i].allele_count > 0) {
			findLongestHaps(v,i);
		}
	}
	return 0;
}


















