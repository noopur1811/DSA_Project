//trie.h
    // Start of include guard
#define TRIE_H

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>

typedef struct TrieNode {
    struct TrieNode *children[4]; // One for each nucleotide: A, C, G, T
    bool isTerminal;        
    struct TrieNode *fail;       
    int key;
} TrieNode;

typedef struct ref_genome {
    int key;
    char *label;
    char *sequence;
    char * significance;
    float frequency;  
} ref_genome;

typedef struct user {
    char *name;
    int age;
    char *gender;
    char *test_sequence;
} user;

ref_genome* hash_t[15]; //declare hash table

// Structure to hold differences and edit distance
typedef struct approx_match{
    int geneKey;           // Reference genome identifier
    int editDistance;      // Calculated edit distance
    int *diffPositions;    // Array of differing positions
    int diffCount;         // Number of differing positions
    int var_pos;
} approx_match;

typedef struct heap{    //heap of structs 'approx_match'
    approx_match *a;           //array of structs
    int size; 
    int r;
}heap;

heap h1;

//function prototypes
TrieNode* create_node();                  
void initialize_trie(TrieNode **root);    // Initializes the Trie with a root node
void insert_into_trie(TrieNode **root, const char *sequence, int key); // Inserts a sequence into the Trie
void free_trie(TrieNode *root);
int getIndex(char nucleotide);
void printtrie_rec(TrieNode *node, const char * reference_genome, int length);
void printtrie(TrieNode * root);
//void print_failure_links(TrieNode *root);
void build_failure_links(TrieNode *root);
int* aho_corasick_search(TrieNode *root, const char *sequence);

//hash func
void report(int *x,const char *s);

//approx_match
int calculate_substitution_distance(const char* query, const char* reference, int* diffPositions);
approx_match find_best_match(const char* query, const char* reference, int key);

//heap
void heap_init(heap *h, int size);
void heap_insert(heap *h, approx_match m);
void report_approx_match(heap *h);
void free_heap(heap *h);


