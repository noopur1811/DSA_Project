#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>

typedef struct variant {
    char gene[30];
    char sequence[150];  
    int chromosome;
    char significance[15];
    float frequency;  
} variant;

variant* hash_t[10] = {NULL};//initialise the hash table to null

void report(int x) {//modify acc to the array of trie return
    variant *v = hash_t[x - 1];  // Adjust for zero-indexed array
    if (v) {
        printf("Gene: %s, Sequence: %s, Significance: %s, Frequency: %.2f\n", v->gene, v->sequence, v->significance, v->frequency);
    } else {
        printf("No variant found at index %d.\n", x);
    }
}

// Tokenize each line from the CSV
void tokenize_line(char* line, variant* rfg) {
    char *p;

    // Parse gene name
    p = strtok(line, ",");
    strncpy(rfg->gene, p, sizeof(rfg->gene));
    rfg->gene[sizeof(rfg->gene) - 1] = '\0';  // Null terminate to avoid overflow

    // Parse sequence
    p = strtok(NULL, ",");
    strncpy(rfg->sequence, p, sizeof(rfg->sequence));
    rfg->sequence[sizeof(rfg->sequence) - 1] = '\0';  // Null terminate to avoid overflow

    // Parse clinical significance
    p = strtok(NULL, ",");
    strncpy(rfg->significance, p, sizeof(rfg->significance));
    rfg->significance[sizeof(rfg->significance) - 1] = '\0';  // Null terminate to avoid overflow

    // Parse frequency of occurrence
    p = strtok(NULL, ",");
    rfg->frequency = atof(p);  // Convert string to float
}

// Read CSV and store variants
int readcsv(const char *filename, variant *rfg, int limit) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Could not open file");
        return 0;
    }

    char line[1024];  // Buffer for each line in the CSV file
    int count = 0;


    // Read each line until EOF or limit is reached
    while (fgets(line, sizeof(line), file) && count < limit) {
        // Remove newline characters
        line[strcspn(line, "\r\n")] = '\0';  // Handle both Windows and Unix-style line endings

        variant *v = &rfg[count];
        tokenize_line(line, v);
        hash_t[count] = v;  // Store in the hash table for quick lookup

        count++;
    }

    fclose(file);
    return count;
}


int main() {
    variant rfg[15];  // Array to hold up to 10 variants
    int num_variants = readcsv("C:/Users/parki/OneDrive/Desktop/DSA/-DSAprojectfilled_reference_genome (2).csv", rfg, 15);
    
    printf("Enter a key from 1-10:\n");
    int x;
    scanf("%d", &x);//input key from the user
    report(x);
    return 0;
}
