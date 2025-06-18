//triee.c

#include "trie.h"
#include <stdbool.h>
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include <limits.h>
#include <math.h>

int getIndex(char nucleotide) {
    if (nucleotide == 'A') return 0;
    if (nucleotide == 'C') return 1;
    if (nucleotide == 'G') return 2;
    if (nucleotide == 'T') return 3;
    return -1; // Invalid character
}


TrieNode* create_node() {
    TrieNode *newNode = (TrieNode *)malloc(sizeof(TrieNode));
    
    // Initialize child pointers to NULL
    for (int i = 0; i < 4; i++) {
        newNode->children[i] = NULL;
    }
    
    // Initialize other properties
    newNode->isTerminal = false;  // Not a terminal node initially
    newNode->key = -1;        // No label initially

    return newNode;
}

void initialize_trie(TrieNode ** root) {
    *root = create_node();   // Create the root node
}
//search full algo time complexity
//O(n * m + q).
//n is number of reference genomes m is length of genomes q is length of query

//insert into trie
void insert_into_trie(TrieNode **root, const char *sequence, int key) {
    // Initialize root node if necessary
    if (*root == NULL) {
        *root = create_node();
    }
    TrieNode *tmp = *root;
    int length = strlen(sequence);

    for (int i = 0; i < length; i++) {
        int index = getIndex(sequence[i]);
        if (index == -1) {
            printf("Invalid character in sequence\n");
            return;
        }
        if (tmp->children[index] == NULL) {
            tmp->children[index] = create_node(); // Create a new node if child is missing
        }
        tmp = tmp->children[index]; // Traverse to the child node
    }

    tmp->isTerminal = true; // Mark the end of a valid sequence

    
    if (tmp->key != -1) {
        tmp->key = -1;  
    } else {
        tmp->key = key;  
    }
}

void build_failure_links(TrieNode *root) {
    TrieNode **queue = malloc(1000 * sizeof(TrieNode *));
  // Increase size for safety, or consider dynamic allocation
    int front = 0, rear = 0;

    // Initialize root's children failure links
    for (int i = 0; i < 4; i++) {
        if (root->children[i]) {
            root->children[i]->fail = root;
            queue[rear++] = root->children[i];
        }
    }

    while (front < rear) {
        TrieNode *current = queue[front++];

        // Set the failure links for each child
        for (int i = 0; i < 4; i++) {
            if (current->children[i]) {
                TrieNode *failure = current->fail;

                // Follow failure links to find a matching child
                while (failure && !failure->children[i]) {
                    failure = failure->fail;
                }

                if (failure) {
                    current->children[i]->fail = failure->children[i];
                    // Inherit terminal properties and key from the failure node if applicable
                    if (failure->children[i]->isTerminal) {
                        current->children[i]->isTerminal = true;
                        current->children[i]->key = failure->children[i]->key;
                    }
                } else {
                    current->children[i]->fail = root;
                }

                queue[rear++] = current->children[i];
            }
        }
    }
    //printf("failure links built\n");
}

// Aho-Corasick search function for DNA sequences
int* aho_corasick_search(TrieNode *root, const char *sequence) {
    int * result = (int *)malloc(2 * sizeof(int));
    TrieNode *current = root;
    int position = 0, x=20, count = 0;
   
    //initialise result with -1
    result[0]=-1;
    result[1]=-1;
    
    while (*sequence) {
        int index = getIndex(*sequence);

        if (index == -1) {
            printf("Invalid nucleotide in query sequence: %c\n", *sequence);
            sequence++;
            position++;
            continue;
        }

        // Follow failure links if no matching child
        while (!current->children[index]) {
            current = current->fail;
        }

        current = current->children[index];

        // Traverse through the output chain if there are matches
        TrieNode *temp = current;
        while (temp != root) {
            if (temp->isTerminal) {
                //printf("Match found for gene %d at position %d\n", temp->key, position);
                result[0] = temp->key;
                result[1] = position;
                count++;
                report(result,sequence);
            }
            temp = temp->fail;
        }

        sequence++;
        position++;
    }
    if(count == 0){
        report(result,sequence);
    }
    return result;
}

void free_trie(TrieNode *root){
    if(root == NULL) return; //trie is empty

    for (int i = 0; i < 4; i++) {
        if (root->children[i] != NULL) {
            free_trie(root->children[i]);
        }
    }

    // Free the label if it exists (for terminal nodes)
    if (root->key != -1) {
        root->key = -1;
    }

    // Finally, free the current node
    free(root);
}

//for debugging 
void printtrie_rec(TrieNode *node, const char * reference_genome, int length){
    // int length = strnlen(reference_genome);
    char new_reference_genome[length+2];
    memcpy(new_reference_genome, reference_genome, length);
    new_reference_genome[length+1] = '\0';

    if(node->isTerminal){
        printf("SEQUENCE: %s, KEY: %d\n", reference_genome, node->key);
    }

    for(int i=0; i<4; i++){
       if (node->children[i] != NULL) {
            // Assign the correct nucleotide character based on the index
            switch (i) {
                case 0: new_reference_genome[length] = 'A'; break;
                case 1: new_reference_genome[length] = 'C'; break;
                case 2: new_reference_genome[length] = 'G'; break;
                case 3: new_reference_genome[length] = 'T'; break;
            }
            // Recursive call with the updated prefix
            printtrie_rec(node->children[i], new_reference_genome, length + 1);
        }
    }
}

void printtrie(TrieNode * root){
    printf("printtrie called\n");
    if (root==NULL){
        printf("TRIE EMPTY\n");
        return;
    }
    char reference_genome[1] = {0};  // Start with an empty reference_genome
    printtrie_rec(root, reference_genome, 0);
}


int calculate_substitution_distance(const char* query, const char* reference, int* diffPositions) {
    int query_len = strlen(query);
    int ref_len = strlen(reference);
    int min_len = query_len < ref_len ? query_len : ref_len;

    int diffCount = 0;

    // Compare characters at aligned positions
    for (int i = 0; i < min_len; i++) {
        if (query[i] != reference[i]) {
            diffPositions[diffCount++] = i;  // Record differing position
        }
    }

    return diffCount;
}

approx_match find_best_match(const char* query, const char* reference, int key) {
    int query_len = strlen(query);  // Length of the query sequence
    int ref_len = strlen(reference); // Length of the reference genome
    int window_size = ref_len;       // The size of the sliding window
    int min_edit_distance = INT_MAX;
    int best_start_pos = -1;

    approx_match best_match = { 
        .editDistance = INT_MAX, 
        .diffPositions = NULL, 
        .diffCount = 0, 
        .var_pos = -1 
    };

    // Sliding window through the query sequence
    for (int i = 0; i <= query_len - window_size; i++) {
        

        // Extract a substring from the query sequence
        char* query_subseq = (char*)malloc((window_size + 1) * sizeof(char));
        strncpy(query_subseq, query + i, window_size);
        query_subseq[window_size] = '\0'; // Null-terminate the substring

       
        // Allocate memory for differing positions
        int* diffPositions = (int*)malloc(window_size * sizeof(int));

        // Use substitution distance to calculate differences
        int substitution_distance = calculate_substitution_distance(query_subseq, reference, diffPositions);
        
        // Check if this match is better (lower substitution distance)
        if (substitution_distance < min_edit_distance) {
            
            min_edit_distance = substitution_distance;
            best_start_pos = i;

            // Free previously stored differing positions if any
            if (best_match.diffPositions != NULL) {
                free(best_match.diffPositions);
            }

            // Update the best match details
            best_match.editDistance = substitution_distance;
            best_match.var_pos = best_start_pos;

            // Count and store differing positions
            int diffCount = 0;
            for (int j = 0; j < window_size && diffPositions[j] != 0; j++) {
                diffCount++;
            }
            best_match.diffCount = diffCount;

            // Allocate memory for best match diffPositions
            best_match.diffPositions = (int*)malloc(diffCount * sizeof(int));
            for (int j = 0; j < diffCount; j++) {
                best_match.diffPositions[j] = diffPositions[j];
               
            }
        }

        // Free temporary memory
        free(query_subseq);
        free(diffPositions);
    }
    best_match.geneKey = key;
    //printf("Debug: Best match found at position %d with edit distance %d\n", best_match.var_pos, best_match.editDistance);
    return best_match;
}



void heap_init(heap *h, int size) {
    h->size = size;
    h->r = -1;
    h->a = (approx_match *)malloc(size * sizeof(approx_match)); // Allocate memory for the heap array
    if (!h->a) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
}


int parent(int i){
    return (i-1)/2;
}
// Insert a new approx_match into the heap
void heap_insert(heap *h, approx_match m) {
    if (h->r == h->size - 1) {                                                          
        printf("Heap is full. Cannot insert.\n");
        return;
    }

    approx_match n; // New match node with given details
    n.geneKey = m.geneKey;
    n.editDistance = m.editDistance;
    n.var_pos = m.var_pos;
    n.diffCount = m.diffCount;
    n.diffPositions = m.diffPositions;

    h->r++; // Insert as rear
    h->a[h->r] = n;
    int i = h->r;

    // Percolate up
    while (i > 0 && h->a[i].editDistance < h->a[(i - 1) / 2].editDistance) {
        approx_match temp = h->a[i];
        h->a[i] = h->a[(i - 1) / 2];
        h->a[(i - 1) / 2] = temp;
        i = (i - 1) / 2;
    }
}

void free_heap(heap *h) {
    if (h->a != NULL) {
        free(h->a);
        h->a = NULL;
        h->r = -1;
        h->size = 0;
    }
}
// Report the closest approx_match from the heap
void report_approx_match(heap *h) {
    if (h->r == -1) {
        printf("The heap is empty.\n");
        return;
    }

    approx_match min = h->a[0];
    ref_genome *rf = hash_t[min.geneKey]; // Retrieve the reference genome with the matching key

    if (rf) {
        printf("Exact matches of any genetic variant were not found, but the genome showcases a risk of:\n");
        printf("The variant of '%s'.\nIt was found to be located at position %d in the input genome.\nThe variant is of the type %s and attacks on chromosome number %d.\nNOTE: there is a high risk of this gene getting carried over to future generations!\n",
               rf->label, min.var_pos, rf->significance, rf->key);
       // printf("%d check edit dist ", min.editDistance);
    } else {
        printf("No valid variant found at the given index %d.\n", min.geneKey);
        return;
    }
   
}

//report
void report(int *x, const char *s){
    //adjust for zero-indexed array
    if (x[0] > 15 ) {
        printf("No variant found at index %d.\n", x[0]);
        return;
    }
    else if(x[0]==-1 && x[1]==-1){
        //printf("No genetic variants detected!\n");
        return;
    }
    ref_genome *v = hash_t[x[0] - 1];
    
    printf("Identified Gene is of the variant name: %s\nSequence of chromosomes in the variant was found to be: %s\nThe gene is of %s type\nFrequency of the identified variant is %.2f\n",
           v->label, v->sequence, v->significance, v->frequency);
    printf("Position at which the genome is affected is %d\n",x[1]-strlen(v->sequence)+2);
    printf("\n");
}
