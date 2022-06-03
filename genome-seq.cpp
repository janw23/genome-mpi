// Author: Jan Wawszczak (jw418479)
// Sequential implementation of genome_index

#include<string>
#include<iostream>
#include<algorithm>
#include<vector>
#include<tuple>
#include<cassert>

#include<random>

int num_occurences_naive(std::string const &genome, std::string const &query) {
    int num_occurences_naive = 0;
    std::size_t pos = 0;
    while ((pos = genome.find(query, pos)) != std::string::npos) {
        num_occurences_naive++;
        pos++;
    }
    return num_occurences_naive;
}

template <typename T>
void printvec(std::string name, std::vector<T> const &vec) {
    return;
    std::cerr << name << ":";
    for (T e : vec) std::cerr << " " << e;
    std::cerr << "\n";
}

// Global-sorting prefix-doubling algorithm.
// Output is a suffix array of genome.
// TODO suffix array type must fit indices of genome. That's why it's size_t here.
std::vector<std::size_t> gspf(std::string const &genome) {
    // Put k-mers into B.
    // TODO Squeeze k>1 - mers into B here to optimize num iterations.
    // TODO starting with k>1 would optimize runtime, but not mem usage.
    std::vector<std::size_t> B(genome.begin(), genome.end());

    std::vector<std::size_t> SA(B.size());

    { // Sort B, keeping original indices.
        std::vector<std::pair<char, std::size_t>> tmp(B.size());
        for (std::size_t i = 0; i < B.size(); i++) {
            tmp[i] = std::make_pair(B[i], i);
        }

        std::sort(tmp.begin(), tmp.end()); // TODO no need to compare by second element of pair
        
        for (std::size_t i = 0; i < tmp.size(); i++) {
            B[i] = tmp[i].first;
            SA[i] = tmp[i].second;
        }
    }

    { // Rebucket
        std::size_t g = 0;
        std::size_t prev = B[0];
        B[0] = 0;
        for (std::size_t i = 1; i < B.size(); i++) {
            if (prev != B[i]) {
                g = i;
            }
            prev = B[i];
            B[i] = g;
        }
    }

    // Iterate algorithm, starting with k chosen for initial k-mers TODO here k=1
    for (std::size_t h = 1;; h *= 2) {
        { // Reorder to string-order. TODO in place?
            auto tmp = B;
            for (std::size_t i = 0; i < B.size(); i++) {
                tmp[SA[i]] = B[i];
            }
            std::swap(tmp, B);
        }

        // Shift B.
        std::vector<std::size_t> B2(B.size());
        std::copy(B.begin() + h, B.end(), B2.begin()); // copy values shifted by h

        {   // Sort according to B, B2 keeping original indices.
            std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> tmp(B.size());
            for (std::size_t i = 0; i < B.size(); i++) {
                tmp[i] = std::make_tuple(B[i], B2[i], i);
            }
            std::sort(tmp.begin(), tmp.end()); // TODO no need to compare by third element of tuple

            for (std::size_t i = 0; i < B.size(); i++) {
                B[i] = std::get<0>(tmp[i]);
                B2[i] = std::get<1>(tmp[i]);
                SA[i] = std::get<2>(tmp[i]);
            }
        }

        // Rebucket.
        {
            std::size_t g = 0;
            std::size_t prev = B[0];
            B[0] = 0;
            std::size_t singletons = 1;
            for (std::size_t i = 1; i < B.size(); i++) {
                if (prev != B[i] || B2[i-1] != B2[i]) { // TODO this seems correct but is different than paper?
                    g = i;
                    singletons++;
                }
                prev = B[i];
                B[i] = g;
            }
            if(singletons == B.size()) break; // Check algorithm-done condition.
        }
    }

    return SA;
}

int num_occurences(std::vector<std::size_t> const &genome_sa, std::string const &genome, std::string const &query) {    
    std::size_t leftmost, rightmost;
    bool leftmost_found = false, rightmost_found = false;
    
    std::size_t l = 0, r = genome.size() - 1, c;
    while (l <= r) {
        c = (l + r) / 2;

        int cmp = genome.compare(genome_sa[c], query.size(), query);
        if (cmp == 0) {
            r = c;
            if (l == r) {
                leftmost = c;
                leftmost_found = true;
                break;
            }
        } else if (cmp > 0) {
            r = c - 1;
        } else {
            l = c + 1;
        }
    }
    if (!leftmost_found) return 0;

    l = leftmost;
    r = genome.size() - 1;
    while (l <= r) {
        c = (l + r + 1) / 2;

        int cmp = genome.compare(genome_sa[c], query.size(), query);
        if (cmp == 0) {
            l = c;
            if (l == r) {
                rightmost = c;
                rightmost_found = true;
                break;
            }
        } else if (cmp > 0) {
            r = c - 1;
        } else {
            l = c + 1;
        }
    }

    assert(rightmost_found);
    return rightmost - leftmost + 1;
}

int main() {
    // std::string genome = "ACGTACACACCCGCTACCGACCGTC$";
    // std::string query = "ACC";  

    // // std::cout << num_occurences_naive(genome, query) << std::endl;

    // auto genome_sa = gspf(genome);
    // assert(genome_sa.size() == genome.size());

    // for (auto i : genome_sa) {
    //     std::cout << i << " ";
    // }
    // std::cout << std::endl;

    // auto num = num_occurences(genome_sa, genome, query);
    // std::cerr << "num occurences: " << num << "\n";
    srand(1);

    for (int test = 0; test < 1000; test++) {
        static const std::string alphabet = "ACGT";
        std::size_t genome_len = 5 + rand() % 10000;
        
        std::string genome;
        genome.reserve(genome_len);

        for (size_t i = 0; i < genome_len; i++) {
            genome.push_back(alphabet[rand() % alphabet.size()]);
        }
        genome.push_back('$');
        auto genome_sa = gspf(genome);

        for (int num_query = 0; num_query < 10; num_query++) {
            std::string query;
            size_t query_len = 1 + rand() % 100;
            query.reserve(query_len);  
            
            for (std::size_t i = 0; i < query_len; i++) {
                query.push_back(alphabet[rand() % alphabet.size()]);
            }

            auto naive = num_occurences_naive(genome, query);
            auto numsa = num_occurences(genome_sa, genome, query);

            assert(naive == numsa);
        }
    }

    return 0;
}