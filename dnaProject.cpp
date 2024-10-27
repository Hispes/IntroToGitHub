#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

bool isValidBase(char base){ // this need some comments

    // The function will determine whether a given character is a valid DNA base. 
    if(base == 'A' || base == 'G' || base == 'C' || base == 'T'){
        return true;
    }
    return false;
}

bool isValidStrand(string strand){
    
    if(strand.empty()){
        return false;
    }
    
    for(size_t i = 0; i < strand.length(); i++){
        if(!isValidBase(strand[i]))
            return false;
    }
        
    return true; 
}


double strandSimilarity(string strand1, string strand2) {
    
    if(strand1.length() != strand2.length()){
        return 0.0;
    }
    
    double matches = 0;
    
    for(size_t i = 0; i < strand1.length(); i++){
        if(strand1[i] ==  strand2[i]){
            matches++;
        }
    }
    
    double similarity = matches/ (strand1.length());
    
    return similarity;
    
}

// int bestStrandMatch(const string& input_strand, const string& target_strand) {
//     const string& longer_strand = (input_strand.length() >= target_strand.length()) ? input_strand : target_strand;
//     const string& shorter_strand = (input_strand.length() < target_strand.length()) ? input_strand : target_strand;

//     double maxScore = 0;
//     int bestInd = 0;
//     int maxIterations = longer_strand.length() - shorter_strand.length();

//     for (int i = 0; i <= maxIterations; i++) {
//         string substring = longer_strand.substr(i, shorter_strand.length());
//         double similarity = strandSimilarity(substring, shorter_strand);

//         if (similarity > maxScore) {
//             maxScore = similarity;
//             bestInd = i;
//         }
//     }
    
//     cout << "Best similarity score: " << maxScore << endl;
//     return bestInd;
// }

// int bestStrandMatch(const string& input_strand, const string& target_strand) {
    
//     const string& longer_strand = (input_strand.length() >= target_strand.length()) ? input_strand : target_strand;
//     const string& shorter_strand = (input_strand.length() < target_strand.length()) ? input_strand : target_strand;
    
    
    
//     double maxScore = 0;
//     int bestInd = 0;
//     int maxIterations = longer_strand.length() - shorter_strand.length();

//     for (int i = 0; i <= maxIterations; i++) {
//         string substring = longer_strand.substr(i, shorter_strand.length());
//         double similarity = strandSimilarity(substring, shorter_strand);

//         if (similarity > maxScore) {
//             maxScore = similarity;
//             bestInd = i;
//         }
//     }
    
//     // Only print the score if it's greater than 0
//     if (maxScore > 0) {
//         cout << "Best similarity score: " << maxScore << endl;
//     } else {
//         cout << "Best similarity score: 0.0" << endl;
//     }
//     return bestInd;
// }



int bestStrandMatch(string input_strand, string target_strand){
    
    if(input_strand.length() < target_strand.length()){
        cout << "Best similarity score: 0.0" << endl;
        return -1;
    }
    
    double maxScore = 0;
    int bestInd = 0;
    int maxIterations = input_strand.length() - target_strand.length();
    
    for(int i = 0; i <= maxIterations; i++){
        
        string substring = input_strand.substr(i, target_strand.length());
        double similarity = strandSimilarity(substring, target_strand);
        
        if(similarity > maxScore){
            maxScore = similarity;
            bestInd = i;
        }
    }
    cout << "Best similarity score: " << maxScore << endl;;
    return bestInd;
}

void identifyMutations(string input_strand, string target_strand) {
    int alignmentIndex = bestStrandMatch(input_strand, target_strand);
    cout << "Best alignment index: " << alignmentIndex << endl;

    // Pad the shorter strand with spaces to align it properly
    if (input_strand.length() < target_strand.length()) {
        input_strand = string(alignmentIndex, ' ') + input_strand + string(target_strand.length() - input_strand.length() - alignmentIndex, ' ');
    } else if (target_strand.length() < input_strand.length()) {
        target_strand = string(alignmentIndex, ' ') + target_strand + string(input_strand.length() - target_strand.length() - alignmentIndex, ' ');
    }

    bool mutation_found = false;
    
    for (size_t i = 0; i < max(input_strand.length(), target_strand.length()); i++) {
        if (i < input_strand.length() && i < target_strand.length()) {
            if (input_strand[i] == ' ' && target_strand[i] != ' ') {
                cout << "Insertion at position " << i + 1 << ": " << target_strand[i] << " is inserted in target strand" << endl;
                mutation_found = true;
            } else if (input_strand[i] != ' ' && target_strand[i] == ' ') {
                cout << "Deletion at position " << i + 1 << ": " << input_strand[i] << " is deleted in target strand" << endl;
                mutation_found = true;
            } else if (input_strand[i] != target_strand[i] && input_strand[i] != ' ' && target_strand[i] != ' ') {
                cout << "Substitution at position " << i + 1 << ": " << input_strand[i] << " -> " << target_strand[i] << endl;
                mutation_found = true;
            }
        }
    }

    if (!mutation_found) {
        cout << "No mutations found." << endl;
    }
}

void transcribeDNAtoRNA(string strand){
    for(unsigned i = 0; i < strand.length(); i++){
        if(strand[i] == 'T'){
            strand[i] = 'U';
        }
    }
    cout << "The transcribed DNA is: " << strand << endl;
}

void reverseComplement(string strand){
    
    
    reverse(strand.begin(), strand.end());
    
    for(unsigned i = 0; i < strand.length(); i++){
        if(strand[i] == 'A'){
            strand[i] = 'T';
        }else if(strand[i] == 'T'){
            strand[i] = 'A';
        }else if(strand[i] == 'C'){
            strand[i] = 'G';
        }else if(strand[i] == 'G'){
            strand[i] = 'C';
        }
    }
    cout << "The reverse complement is: " << strand << endl;
}


void getCodingFrames(string strand) {
    bool foundFrame = false; 

    
    for (size_t i = 0; i < strand.length(); i++) {
        
        if (strand[i] == 'A' && strand[i + 1] == 'T' && strand[i + 2] == 'G') {
            
            
            for (size_t j = i + 3; j < strand.length(); j += 3) { // Move in steps of 3 letters
                
                string codon = strand.substr(j, 3); 

                
                if (codon == "TAA" || codon == "TAG" || codon == "TGA") {
                    // Print the whole sequence from "ATG" to the stop codon
                    for (size_t k = i; k <= j + 2; k++) {
                        cout << strand[k];
                    }
                    cout << endl;
                    foundFrame = true;
                    
                    
                    i = j + 2; 
                    break; 
                }
            }
        }
    }

   
    if (!foundFrame) {
        cout << "No reading frames found." << endl;
    }
}



int main(){
    string dna1, dna2;
    double score = 0;
    int userInput;
    bool runProgram = true;
    while(runProgram){
        cout << "--- DNA Analysis Menu ---" << endl;
        cout << "1. Calculate the similarity between two sequences of the same length\n2. Calculate the best similarity between two sequences of either equal or unequal length\n3. Identify mutations\n4. Transcribe DNA to RNA\n5. Find the reverse complement of a DNA sequence \n6. Extract coding frames \n7. Exit" << endl;
        cout << "Please enter your choice (1 - 7):" << endl;
        cin >> userInput;

        switch(userInput){
            case 1:
                cout << "Enter the first DNA sequence: " << endl;
                cin >> dna1;
                if(!isValidStrand(dna1)){
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the first DNA sequence: " << endl;
                    cin >> dna1;
                }
                cout << "Enter the second DNA sequence:" << endl;
                cin >> dna2;
                if(!isValidStrand(dna2)){
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the second DNA sequence:" << endl;
                    cin >> dna2;
                }
                score = strandSimilarity(dna1,dna2);
                if(score == 0){
                    cout << "Error: Input strands must be of the same length." << endl;
                }else{
                    cout << "Similarity score: " <<  score << endl;
                }
                break;
            case 2:
                cout << "Enter the first DNA sequence: " << endl;
                cin >> dna1;
                if(!isValidStrand(dna1)){
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the first DNA sequence: " << endl;
                    cin >> dna1;
                }
                cout << "Enter the second DNA sequence:" << endl;
                cin >> dna2;
                if(!isValidStrand(dna2)){
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the second DNA sequence:" << endl;
                    cin >> dna2;
                }
                bestStrandMatch(dna1,dna2);
                break;
            case 3:
                cout << "Enter the first DNA sequence: " << endl;
                cin >> dna1;
                cout << "Enter the second DNA sequence:" << endl;
                cin >> dna2;
                identifyMutations(dna1,dna2);
                break;
            case 4:
                cout << "Enter the DNA sequence to be transcribed: " << endl;
                cin >> dna1;
                if(!isValidStrand(dna1)){
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the DNA sequence to be transcribed: " << endl;
                    cin >> dna1;
                }
                transcribeDNAtoRNA(dna1);
                break;
            case 5:
                cout << "Enter the DNA sequence: " << endl;
                cin >> dna1;
                if(!isValidStrand(dna1)){
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the DNA sequence: " << endl;
                    cin >> dna1;
                }
                reverseComplement(dna1);
                break;
            case 6:
                cout << "Enter the DNA sequence: " << endl;
                cin >> dna1;
                if(!isValidStrand(dna1)){
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the DNA sequence: " << endl;
                    cin >> dna1;
                }
                cout << "The extracted reading frames are: " << endl;
                getCodingFrames(dna1);
                break;
            case 7:
                cout << "Exiting program." << endl;
                runProgram = false;
                break;
            default:
                cout << "Invalid input. Please select a valid option." << endl;
                break;
                cout << "--- DNA Analysis Menu ---" << endl;
                cout << "1. Calculate the similarity between two sequences of the same length\n2. Calculate the best similarity between two sequences of either equal or unequal length\n3. Identify mutations\n4. Transcribe DNA to RNA\n5. Find the reverse complement of a DNA sequence \n6. Extract coding frames \n7. Exit" << endl;
                cout << "Please enter your choice (1 - 7):" << endl;
                cin >> userInput;
        }
    }
    return 0;
}