#include "trie.h"
#include <stdbool.h>
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include <fcntl.h>
#include <unistd.h>


#define LINESIZE 1024
//csv 1 functions 
// Reads a single line from a file descriptor (fd)
int readline(int fd, char *line, int linesize) {
    char ch, *p = line;
    int count = 0;
    while (count < linesize && read(fd, &ch, 1)) {
        if (ch != '\n')  {
            *line++ = ch;
            count++;
        } else {
            break;
        }
    }
    *line = 0;
    return line - p;
}

// Removes commas from a CSV string
void remove_comma(char *p) {
    char *q;
    while (*p) {
        if (*p == ',') {
            q = p;
            while (*q) {
                *q = *(q + 1);
                q++;
            }    
        } else {
            p++;
        }
    }
}

// Tokenizes a line into a `ref_genome` struct
void tokenize_line(char* line, ref_genome* rfg) {
    char *p = strtok(line, ",");  
    if (p != NULL) {
        rfg->key = atoi(p);  // Convert the first token to `key`
    }

    p = strtok(NULL, ",");  // Get the second token (label)
    if (p != NULL) {
        rfg->label = strdup(p);  // Copy the label
    }

    p = strtok(NULL, ",");  // Get the third token (sequence)
    if (p != NULL) {
        rfg->sequence = strdup(p);  // Copy the sequence
    }

    p = strtok(NULL, ",");  // Get the fourth token (significance)
    if (p != NULL) {
        rfg->significance = strdup(p);
    }

    p = strtok(NULL, ",");  // Get the fifth token (frequency)
    if (p != NULL) {
        rfg->frequency = atof(p);
    }
}

int readcsv(const char *filename, ref_genome *rfg, int limit) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Could not open file");
        return 0;
    }

    char line[LINESIZE];
    int count = 0;

    while (fgets(line, sizeof(line), file) && count < limit) {
        // Remove trailing newline if it exists
        line[strcspn(line, "\r\n")] = 0;

        // Check for basic CSV structure
        if (strchr(line, ',') == NULL) {
            printf("Skipping malformed line: %s\n", line);
            continue;
        }

        // Tokenize and store in rfg[count]
        tokenize_line(line, &rfg[count]);
        ref_genome *v =&rfg[count]; //add the values in hash table consecutively
        hash_t[count] = v;
        count++;
    }

    fclose(file);
    return count;
}


// Function to print the contents of ref_genome structs (for debugging purposes)
void print_ref_genomes(ref_genome *rfg, int count) {
    printf("Printing reference genomes:\n");
    for (int i = 0; i < count; i++) {
        printf("Key: %d\n", rfg[i].key);
        printf("Label: %s\n", rfg[i].label);
        printf("Sequence: %s\n", rfg[i].sequence);
        printf("Significance: %s\n", rfg[i].significance);
        printf("Frequency: %.2f\n\n", rfg[i].frequency); // Assuming frequency has 2 decimal places
    }
}


void insert_all_into_trie(TrieNode *root, ref_genome *rfg, int count) {
    for (int i = 0; i < count; i++) {
        // Insert the sequence and label into the Trie
        insert_into_trie(&root, rfg[i].sequence, rfg[i].key);
    }
    //printf("sequences inserted in trie\n");
}

//functions for second csv
int load_test_cases(const char *filename, user *cases, int limit) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Could not open file");
        return 0;
    }
    
    int count = 0;
    char line[350];
    while (fgets(line, sizeof(line), file) && count < limit) {
        user *test = &cases[count];

        // Allocate memory for each string
        test->name = malloc(50 * sizeof(char)); // Adjust size as needed
        test->gender = malloc(10 * sizeof(char)); // Adjust size as needed
        test->test_sequence = malloc(200 * sizeof(char)); // Adjust size as needed

        if (!test->name || !test->gender || !test->test_sequence) {
            perror("Failed to allocate memory");
            fclose(file);
            return count; // Return the number of successfully loaded cases
        }

        sscanf(line, "%49[^,],%d,%9[^,],%149[^\n]", test->name, &test->age, test->gender, test->test_sequence);
        count++;
    }
    
    fclose(file);
    return count;
}

void free_test_cases(user *cases, int count) {
    for (int i = 0; i < count; i++) {
        free(cases[i].name);
        free(cases[i].gender);
        free(cases[i].test_sequence);
    }
}

//Find edit distance- When there is no exact match
//report func



int main(){
    TrieNode *root = NULL; 
    user cases[35];
    heap_init(&h1,150);
    int num_cases = load_test_cases("trial.csv", cases, 20);
    int *match_array = (int*)malloc(2*sizeof(int));

    initialize_trie(&root);  //using csv file

    ref_genome rfg[100];  // Assuming a maximum of 100 records
    int limit = 100;  // Limit of records to read

    
    int count = readcsv("filled_reference_genome (2).csv", rfg, limit); // Read the CSV file and populate the array
    //printf("Records read: %d\n", count);  // Debug print to check if we read records

    insert_all_into_trie(root, rfg, count);
    //printtrie(root);

    build_failure_links(root);
    //printtrie(root);


    //second csv trial 
    if (num_cases == 0) {
        printf("No test cases loaded.\n");
        return 1;
    }

    for (int i = 0; i < num_cases; i++) {
        printf("USER: %d\n Name: %s\n Age: %d\n Gender: %s\n",i+1,cases[i].name,cases[i].age,cases[i].gender);
        printf("\nReport:\n");

        match_array = aho_corasick_search(root, cases[i].test_sequence);

        if(match_array[0] == -1 && match_array[1] == -1){
       

            for(int j = 0; j<count; j++){
            
                approx_match app_match= find_best_match(cases[i].test_sequence, rfg[j].sequence,j);
                heap_insert(&h1,app_match);
            }

            report_approx_match(&h1);
            free_heap(&h1);
            heap_init(&h1,15);
        }
        printf("\n");
        puts("---------------------------------------------------------------------------------------------");
    }

   
    free_test_cases(cases, num_cases);          // Free the allocated memory

    for (int i = 0; i < count; i++) {           // Free memory allocated for labels and sequences
        free(rfg[i].label);
        free(rfg[i].sequence);
        free(rfg[i].significance);
    }
}