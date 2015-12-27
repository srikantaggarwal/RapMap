#include <sdsl/suffix_arrays.hpp>
#include <fstream>
#include <ostream>
#include <iostream> 
#include <map>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include <sys/time.h>
#include <ctype.h>

#define K 8

using namespace std;
using namespace sdsl;

// Used to measure intervals and absolute times
typedef int64_t msec_t;

// Get current time in milliseconds from the Epoch (Unix)
// or the time the system started (Windows).
msec_t time_ms(void);

class RapMapCSA {

private:
	const char *concatTxtPtr;
	csa_wt<> csa;
	uint64_t csaSize;
	uint64_t stringSize;
	const char EOL = 'A'-1;

public:
	RapMapCSA(const char *str,string outputDir) {
		concatTxtPtr = str;
	        construct_config::byte_algo_sa = sdsl::SE_SAIS;
		string filename = outputDir + "csatemp.txt";
		ofstream output(filename);
		output << concatTxtPtr;
		output.close();
   	        construct(csa, filename, 1);
		csaSize = csa.size();
		stringSize = csaSize-1;
	}

	void populateHashMap(map<uint16_t, pair<uint32_t,uint32_t> > &kmerMap)
	{
     	  	uint64_t it=0;
		uint16_t num = 0;
		uint64_t cur = 0;	 

		while (cur < stringSize)  {        // loop getting single characters
			num = num << 2;
			num = num & ~(0xFFFF<<(2*K)); 
			switch(concatTxtPtr[cur])
			{
				case 'a':
				case 'A':break;
				case 'C':
				case 'c': num += 1;
					  break;	
				case 'g':
				case 'G': num += 2;
					  break;
				case 'T':
				case 't': num += 3;
					  break;	
				default:  break; // treated as an A
						
			}
			cur++;
			if(cur >= K)
			{
				if(kmerMap.find(num) == kmerMap.end())
				{
					uint64_t left{0};
					uint64_t right{0};
					backward_search(csa,0,csaSize-1,concatTxtPtr+it,concatTxtPtr+it+K,left,right);
					kmerMap[num] = make_pair((uint32_t)left,(uint32_t)right);
				}
				it++;	
			}
		}
		
	}

	uint32_t getPositionAt(uint32_t index)
	{
		return csa[index];
	}

	pair<char, char> getCharsPSI(uint32_t index1, uint32_t index2, const char *read, uint32_t readSize) {
		char ch1 = EOL;
		char ch2 = EOL;
		
		if (index1 != 0) {
			for (int i = 0; i < csa.C.size(); ++i) {
				int start = csa.C[i];
				int end = csa.C[i+1];
				if (start <= index1 && index1 < end) {
					ch1 = csa.comp2char[i];
					break;
				}
			}
		}
		if (index2 != readSize) {
			ch2 = read[index2];
		}				
	
		return make_pair(ch1, ch2);
	}
		
	pair<uint64_t, uint64_t> getMMPIntervalPSI(uint64_t left,uint64_t right,uint32_t startAt,const char *read,uint32_t readSize)
	{
		uint32_t left_=left;
		uint32_t right_=right;
		uint32_t maxMMP{0};

		while(left_ <= right_)
		{
			uint32_t middle = (left_+right_)/2;
			uint32_t curMMP = 0;

			uint32_t index = middle;			
			pair<char, char> chars = getCharsPSI(index, startAt+curMMP, read, readSize);
				
			while (chars.second != EOL && chars.first == chars.second) {
				curMMP++;
				index = csa.psi[index];
				chars = getCharsPSI(index, startAt+curMMP, read, readSize);
			}

			if (curMMP > maxMMP)
				maxMMP = curMMP;

			if (chars.second == EOL)
				break;
							
			if (chars.first < chars.second) {
				left_ = middle+1;
			}
			else {
				right_ = middle-1;
			}
		}
				
		uint64_t l{0};
		uint64_t r{0};
		string mmp;
		for (int i = startAt; i < startAt+maxMMP; ++i) {
			char ch = read[i];
			mmp += ch;
		}
		cout << "MMP is " << mmp << endl;	
                backward_search(csa, 0, csaSize-1, mmp.begin(), mmp.end(), l, r);
                return make_pair(l, r);
	}

	pair<char, char> getCharsNaive(uint32_t index1, uint32_t index2, const char *read, uint32_t readSize) {
		char ch1 = EOL;
		char ch2 = EOL;

		if (index1 < stringSize) {
			ch1 = concatTxtPtr[index1];
		}
		if (index2 < readSize) {
			ch2 = read[index2];
		}

		return make_pair(ch1, ch2);
	}

	pair<uint64_t, uint64_t> getMMPIntervalNaive(uint64_t left,uint64_t right,uint32_t startAt,const char *read,uint32_t readSize)
        {
                uint32_t maxMMP{0};
        
	        uint32_t left_ = left; 
	        uint32_t leftPos = csa[left_];
        
		uint32_t LP = K;
		pair<char, char> chars = getCharsNaive(leftPos+LP, startAt+LP, read, readSize);
		while (chars.second != EOL && chars.first == chars.second) {
			LP++;
			chars = getCharsNaive(leftPos+LP, startAt+LP, read, readSize);
		}

		if (maxMMP < LP) {
			maxMMP = LP;
		}

		if (chars.second != EOL) {	
                	uint32_t right_=right;
                	uint32_t rightPos = csa[right_];
			
			uint32_t RP = K;
			chars = getCharsNaive(rightPos+RP, startAt+RP, read, readSize);
			while (chars.second != EOL && chars.first == chars.second) {
				RP++;
				chars = getCharsNaive(rightPos+RP, startAt+RP, read, readSize);
			}

			if (maxMMP < RP) {
				maxMMP = RP;
			}

			if (chars.second != EOL) {
                		while(left_ <= right_)
				{
					uint32_t middle = (left_+right_)/2;
					uint32_t midPos = csa[middle];
					uint32_t MP = min(LP, RP);
					
					chars = getCharsNaive(midPos+MP, startAt+MP, read, readSize);
					while (chars.second != EOL && chars.first == chars.second) {
						MP++;
						chars = getCharsNaive(midPos+MP, startAt+MP, read, readSize);
					}

					if (maxMMP < MP) {
						maxMMP = MP;
					}

					if (chars.second == EOL) {
						break;
					} 
					
					if (chars.first < chars.second) {
						left_ = middle+1;
						LP = MP;
					}
					else {
						right_ = middle-1;
						RP = MP;
					}
				 }
			}	
		}
			
                uint64_t l{0};
                uint64_t r{0};
		string mmp;
		for (int i = startAt; i < startAt+maxMMP; ++i) {
			char ch = read[i];
			mmp += ch;
		}
		cout << "MMP is " << mmp << endl;	
                backward_search(csa, 0, csaSize-1, mmp.begin(), mmp.end(), l, r);
                return make_pair(l, r);
        }
};

msec_t time_ms(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (msec_t)tv.tv_sec * 1000 + tv.tv_usec / 1000;
}

int getHashKey(string str) 
{
	int size = (int)str.size();
	int i = 0;
	int num = 0;
	while (i < size) {
		num = num<<2;
		switch(str[i])
		{
			case 'a':
			case 'A':break;
			case 'C':
			case 'c': num += 1;
				  break;	
			case 'g':
			case 'G': num += 2;
				  break;
			case 'T':
			case 't': num += 3;
				  break;	
			default:  break; // treated as an A
					
		}
		i++;
	}
	return num;
}

int main() {
	
	ifstream stringfile ("concatenated_string.txt");
	if (stringfile.is_open())
	{
		string line;
		while ( getline (stringfile,line) )
		{
			transform(line.begin(), line.end(), line.begin(), ::toupper);
			string outputDir = "./";
  			RapMapCSA rmc(line.c_str(),outputDir);
			
			map<uint16_t,pair<uint32_t,uint32_t> > hashMap;
			rmc.populateHashMap(hashMap);
			
			ifstream readsfile ("r1.fq");
			msec_t startTime = time_ms();
			if (readsfile.is_open()) {
				string temp;
				uint64_t count = 0;
				while ( getline(readsfile, temp) ) {
					count %= 4;
					if (count != 1) {
						count++; 
						continue;
					}
					count++;
					string read;
					for (int i = 0; i < temp.size(); ++i) {
						switch(temp[i]) {
							case 'A':
							case 'a':
								read += 'A';
								break;
							case 'C':
							case 'c':
								read += 'C';
								break;
							case 'G':
							case 'g':
								read += 'G';
								break;
							case 'T':
							case 't':
								read += 'T';
								break;
						}
					}	
					cout << read << endl; 
					
					string kMer = read;
					if (read.size() >= 8) {
						kMer = read.substr(0, 8);
					}

					int key = getHashKey(kMer);
					cout << "8-mer is " << kMer << endl;
			
					uint32_t startAt = 0;

					if (hashMap.find(key) == hashMap.end())
						continue;
	
					pair<uint32_t, uint32_t> interval1 = hashMap[key];
					cout << "Initial SA interval is (" << interval1.first << ", " << interval1.second << ")" << endl;
					//pair<uint64_t, uint64_t> interval2 = rmc.getMMPIntervalPSI((uint64_t)interval1.first,(uint64_t)interval1.second,(uint32_t)startAt,read.c_str(),read.size());
					pair<uint64_t, uint64_t> interval2 = rmc.getMMPIntervalNaive((uint64_t)interval1.first,(uint64_t)interval1.second,(uint32_t)startAt,read.c_str(),read.size());
					cout << "(" << interval2.first << ", " << interval2.second << ")" << endl;				
	
					int pos1 = rmc.getPositionAt(interval2.first);
					int pos2 = rmc.getPositionAt(interval2.second);
					//Need to add check to ensure NIP doesn't exceed read length
					int NIP = startAt+K;
					while(NIP < line.size()) {
						char ch1 = line[pos1+NIP];
						char ch2 = line[pos2+NIP];
						if (ch1 != ch2)
							break;
						NIP++;
					}	
					cout << "NIP occurs at position " << NIP << endl;
					
					cout << "Left end of interval is " << line.substr(pos1, NIP+5-startAt+1) << endl;
					cout << "Right end of interval is " << line.substr(pos2, NIP+5-startAt+1) << endl;
				}
				readsfile.close();
			}
			msec_t endTime = time_ms();
			cout << startTime << ", " << endTime << ", " << endTime-startTime << endl;
		}
		
		stringfile.close();
	}

	return 0;
}
