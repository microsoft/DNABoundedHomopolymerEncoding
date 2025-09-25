/* C++ program to convert a binary string to a string over alphabet 4 with bounded homopolymer length.

Install GMP library on linux with: sudo apt-get install libgmp3-dev
For Windows see: https://gmplib.org/
https://cs.nyu.edu/~exact/core/gmp/index.html also has some instructions for Windows but might be outdated. 
To run a C++ program with GMP library on linux use: 
g++ -O3 mycppprog.cpp -lgmpxx -lgmp -o mycppprog

**************************************************************************/
#include <iostream>
#include <chrono>
#include <string.h>
#include <vector>
#include <gmpxx.h>
using namespace std;
using namespace std::chrono;

class BoundedHomopolymerEncoder{
public:
	BoundedHomopolymerEncoder(int x,int y,int z):max_homopolymer_run_length(x),encoding_length(y),input_data_length(z){

		initialize_FSM(max_homopolymer_run_length);
		initialize_number_paths(encoding_length);

		max_data_len=number_paths[encoding_length][0].get_str(2).length()-1;
		cout<<"Max data bits that can be encoded: "<<max_data_len<<endl;

		if(input_data_length>max_data_len){
			cout<<"ERROR: input_data is too long!"<<endl;
			return;
		}
	}

	string encode(string input_data){
		if (max_homopolymer_run_length==1) return encode_nohomopolymer(input_data);
		return find_Nth_string(lint(input_data,2));
	}

	int max_data_length(){
		return max_data_len;
	}
private:

	typedef mpz_class lint; //long integer
	typedef vector<vector<lint>> matrix_lint;


	int number_states;
	int size_alphabet;
	int encoding_length;
	int max_homopolymer_run_length;
	int input_data_length;
	int max_data_len;

	const vector<vector<int>> FSM1{{1,2,3,4},{-1,2,3,4},{1,-1,3,4},{1,2,-1,4},{1,2,3,-1}};

	const vector<vector<int>> FSM2{{1,2,3,4},{5,2,3,4},{1,6,3,4},{1,2,7,4},{1,2,3,8},{-1,2,3,4},{1,-1,3,4},{1,2,-1,4},{1,2,3,-1}};

	const vector<vector<int>> FSM3{{1,2,3,4},{5,2,3,4},{1,6,3,4},{1,2,7,4},{1,2,3,8},{9,2,3,4},{1,10,3,4},{1,2,11,4},{1,2,3,12},{-1,2,3,4},{1,-1,3,4},{1,2,-1,4},{1,2,3,-1}};

	const vector<vector<int>> FSM4{{1,2,3,4},{5,2,3,4},{1,6,3,4},{1,2,7,4},{1,2,3,8},{9,2,3,4},{1,10,3,4},{1,2,11,4},{1,2,3,12},{13,2,3,4},{1,14,3,4},{1,2,15,4},{1,2,3,16},{-1,2,3,4},{1,-1,3,4},{1,2,-1,4},{1,2,3,-1}};

	const vector<vector<int>> FSM5{{1,2,3,4},{5,2,3,4},{1,6,3,4},{1,2,7,4},{1,2,3,8},{9,2,3,4},{1,10,3,4},{1,2,11,4},{1,2,3,12},{13,2,3,4},{1,14,3,4},{1,2,15,4},{1,2,3,16},{17,2,3,4},{1,18,3,4},{1,2,19,4},{1,2,3,20},{-1,2,3,4},{1,-1,3,4},{1,2,-1,4},{1,2,3,-1}};


	vector<vector<int>> FSM; //Rows correspond to states, columns correspond to alphabet and the entries correspond to next state. FSM[i][a] is the next state if current state is 'i' and next symbol is 'a'. It is set to -1 if the transition is not allowed.

	matrix_lint number_paths; //number_paths(t,s)= number of t-length paths starting from state s


	void initialize_FSM(int max_homopolymer_run){
		switch (max_homopolymer_run)
		{
			case 1: FSM=FSM1; break;
			case 2: FSM=FSM2; break;
			case 3: FSM=FSM3; break;
			case 4: FSM=FSM4; break;
			case 5: FSM=FSM5; break;
			default: cout<<"ERROR: Choose max homopolymer run length between 1 and 5"<<endl;
		}
		number_states=FSM.size();
		size_alphabet=FSM[0].size();
		
		return;
	}

	void initialize_number_paths(int T){
		number_paths=matrix_lint(T+1,vector<lint>(number_states,1)); //default set to 1
		for(int t=1;t<=T;t++){
			for(int s=0;s<number_states;s++){
				lint count=0;
				for(int sigma=0;sigma<size_alphabet;sigma++){
					int next_state=FSM[s][sigma];
					if (next_state!=-1) count+=number_paths[t-1][next_state];
				}
				number_paths[t][s]=count;
			}
		}
	}


	string find_Nth_string(lint N){
	 	// Finds the Nth path in the list of all paths of length encoding_length ordered in lexicographic order.
		// Indexing starts from 0.
		string answer;
		answer.resize(encoding_length);
		// if (number_paths[encoding_length][0]<=N) {
		// 	cout<<"ERROR in encoding"<<endl;
		// 	return answer;
		// }

		int current_state=0; // Begin at the start state
		

		lint count=N; // count = N - number of paths strictly below current partial answer
		for (int t=1;t<=encoding_length;t++){
			for(int sigma=0;sigma<size_alphabet;sigma++){
				int next_state=FSM[current_state][sigma];
				if (next_state==-1) continue;
				if (count<number_paths[encoding_length-t][next_state]){
					answer[t-1]='0'+sigma;
					current_state=next_state;
					break;
				}
				else{
					count-=number_paths[encoding_length-t][next_state];
				}
			}
		}
		return answer;
	}

	string encode_nohomopolymer(string input_data){
		//input_data is a string in base 2
		char first_base=2*(input_data[0]-'0')+input_data[1];


		string input_base3=pad_zeros(lint(input_data.substr(2),2).get_str(3),encoding_length-1);
		string encoding;
		encoding.resize(encoding_length);
		encoding[0]=first_base;
		int a=first_base;
		for(int i=1;i<encoding_length;i++){
			a=(a+input_base3[i-1]-'0'+1)%4;
			encoding[i]=a+'0';
		}
		return encoding;
	}

	string pad_zeros(string s,int L){
		//pad zeros at the beginning of string s until it has length L. Assumes that the length is at most L.
		if (s.length()==L) {return s;}
		return string(L-s.length(),'0')+s;
	}
};



class BoundedHomopolymerDecoder{
public:
	BoundedHomopolymerDecoder(int x,int y,int z):max_homopolymer_run_length(x),encoding_length(y),input_data_length(z){
		initialize_FSM(max_homopolymer_run_length);
		initialize_number_paths(encoding_length);
	}

	string decode(string encoded_data){
		if (max_homopolymer_run_length==1) return decode_nohomopolymer(encoded_data);
		return pad_zeros(find_string_position(encoded_data).get_str(2),input_data_length);
	}

	

private:

	typedef mpz_class lint; //long integer
	typedef vector<vector<lint>> matrix_lint;


	int number_states;
	int size_alphabet;
	int encoding_length;
	int input_data_length;
	int max_homopolymer_run_length;

	const vector<vector<int>> FSM1{{1,2,3,4},{-1,2,3,4},{1,-1,3,4},{1,2,-1,4},{1,2,3,-1}};

	const vector<vector<int>> FSM2{{1,2,3,4},{5,2,3,4},{1,6,3,4},{1,2,7,4},{1,2,3,8},{-1,2,3,4},{1,-1,3,4},{1,2,-1,4},{1,2,3,-1}};

	const vector<vector<int>> FSM3{{1,2,3,4},{5,2,3,4},{1,6,3,4},{1,2,7,4},{1,2,3,8},{9,2,3,4},{1,10,3,4},{1,2,11,4},{1,2,3,12},{-1,2,3,4},{1,-1,3,4},{1,2,-1,4},{1,2,3,-1}};

	const vector<vector<int>> FSM4{{1,2,3,4},{5,2,3,4},{1,6,3,4},{1,2,7,4},{1,2,3,8},{9,2,3,4},{1,10,3,4},{1,2,11,4},{1,2,3,12},{13,2,3,4},{1,14,3,4},{1,2,15,4},{1,2,3,16},{-1,2,3,4},{1,-1,3,4},{1,2,-1,4},{1,2,3,-1}};

	const vector<vector<int>> FSM5{{1,2,3,4},{5,2,3,4},{1,6,3,4},{1,2,7,4},{1,2,3,8},{9,2,3,4},{1,10,3,4},{1,2,11,4},{1,2,3,12},{13,2,3,4},{1,14,3,4},{1,2,15,4},{1,2,3,16},{17,2,3,4},{1,18,3,4},{1,2,19,4},{1,2,3,20},{-1,2,3,4},{1,-1,3,4},{1,2,-1,4},{1,2,3,-1}};


	vector<vector<int>> FSM; //Rows correspond to states, columns correspond to alphabet and the entries correspond to next state. FSM[i][a] is the next state if current state is 'i' and next symbol is 'a'. It is set to -1 if the transition is not allowed.

	matrix_lint number_paths; //number_paths(t,s)= number of t-length paths starting from state s


	void initialize_FSM(int max_homopolymer_run){
		switch (max_homopolymer_run)
		{
			case 1: FSM=FSM1; break;
			case 2: FSM=FSM2; break;
			case 3: FSM=FSM3; break;
			case 4: FSM=FSM4; break;
			case 5: FSM=FSM5; break;
			default: cout<<"ERROR: Choose max homopolymer run length between 1 and 5"<<endl;
		}
		number_states=FSM.size();
		size_alphabet=FSM[0].size();
		
		return;
	}

	void initialize_number_paths(int T){
		number_paths=matrix_lint(T+1,vector<lint>(number_states,1)); //default set to 1
		for(int t=1;t<=T;t++){
			for(int s=0;s<number_states;s++){
				lint count=0;
				for(int sigma=0;sigma<size_alphabet;sigma++){
					int next_state=FSM[s][sigma];
					if (next_state!=-1) count+=number_paths[t-1][next_state];
				}
				number_paths[t][s]=count;
			}
		}
	}


	lint find_string_position(string s){
		int T=s.length();
		lint count=0; 
		int current_state=0;
		for(int t=1;t<=T;t++){
			int a = s[t-1]-'0';
			for(int sigma=0;sigma<a;sigma++){
				int next_state=FSM[current_state][sigma];
				if (next_state==-1) continue;
				count+=number_paths[T-t][next_state];
			}
			current_state=FSM[current_state][a];
		}
		return count;

	}

	string decode_nohomopolymer(string encoding){
		//input_data is a string in base 2
		string bit1bit2;
		switch(encoding[0]){
			case '0': bit1bit2="00";break;
			case '1': bit1bit2="01";break;
			case '2': bit1bit2="10";break;
			case '3': bit1bit2="11";break;
			default: cout<<"ERROR in encoding"<<endl;
		}

		string input_base3;
		input_base3.resize(encoding_length-1);
	
		for(int i=0;i<encoding_length-1;i++){
			int shift=encoding[i+1]-encoding[i];
			if(shift<0) shift+=4;
			input_base3[i]=shift-1+'0';
		}
		return bit1bit2+pad_zeros(lint(input_base3,3).get_str(2),input_data_length-2);
	}


	string pad_zeros(string s,int L){
		//pad zeros at the beginning of string s until it has length L. Assumes that the length is at most L.
		if (s.length()==L) {return s;}
		return string(L-s.length(),'0')+s;
	}
};

int main(int argc, char* argv[]){
	srand (time(NULL));
	

	if(argc!=5){
		cout<<"Give max_homopolymer_run_length, encoding_length, input_data_length, number_trials as arguments. input_data should be in binary."<<endl;
		return 0;
	}
	
	int max_homopolymer_run_length=atoi(argv[1]);
	int encoding_length=atoi(argv[2]);
	int input_data_length=atoi(argv[3]);
	int number_trials=atoi(argv[4]);

	BoundedHomopolymerEncoder BHE(max_homopolymer_run_length,encoding_length,input_data_length);

	if(BHE.max_data_length()<input_data_length){
		cout<<"ERROR: Input data too long. Max data length is: "<<BHE.max_data_length()<<endl; 
		return 1;
	}

	BoundedHomopolymerDecoder BHD(max_homopolymer_run_length,encoding_length,input_data_length);

	string input_data[number_trials];

	for (int i=0;i<number_trials;i++){
		input_data[i].resize(input_data_length);
		for (int j=0;j<input_data_length;j++){
			input_data[i][j]='0'+rand()%2;
		}
	}


	string encoded_data[number_trials];

	auto start = high_resolution_clock::now(); 
	
	for (int i=0;i<number_trials;i++){
		encoded_data[i]=BHE.encode(input_data[i]);

		//cout<<"Encoding: "<<encoded_data[i]<<endl;
	}
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start); 
    cout << "Encoding time: "<< duration.count() << " milliseconds" << endl; 


	start = high_resolution_clock::now(); 
	
	for (int i=0;i<number_trials;i++){
		string decoded_data=BHD.decode(encoded_data[i]);
		if (decoded_data!=input_data[i]) cout<<"ERROR: Decoding failed!\n"<<input_data[i]<<endl;
		//cout<<"Decoding: "<<decoded_data<<endl;
	}
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start); 
    cout << "Decoding time: "<< duration.count() << " milliseconds" << endl; 

}
