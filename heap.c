#include<stdio.h>
#include<stdlib.h>

typedef struct variant {
    char gene[30];
    char sequence[150];  
    int chromosome;
    char significance[15];
    float frequency;  
} variant;

variant* hash_t[10] = {NULL};

//new funcs start here
typedef struct match{
    int pos;            //position of longest substring
    int var;          //variant it is closest to
    int ed;             //edit distance
}match;

typedef struct heap{    //heap of structs 'match'
    match *a;           //array of structs
    int size; 
    int r;
}heap;

void heap_init(heap *h, int size) {
    h->size = size;
    h->r = -1;
    h->a = (match *)malloc(size * sizeof(match)); // Allocate memory for the heap array
    if (!h->a) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
}


int parent(int i){
    return (i-1)/2;
}

void heap_insert(heap *h, int variant, int edit, int posn) {
    if (h->r == h->size - 1) {
        printf("Heap is full. Cannot insert.\n");
        return;
    }

    match n; // new match node with given details
    n.ed = edit;
    n.pos = posn;
    n.var = variant;

    h->r++; // Inserted as rear
    h->a[h->r] = n;
    int i = h->r;

    
    while (i > 0 && h->a[i].ed < h->a[parent(i)].ed) {
        match temp = h->a[i];
        h->a[i] = h->a[parent(i)];
        h->a[parent(i)] = temp;
        i = parent(i); 
    }
}


void report_2(heap *h) {
    if (h->r == -1) {
        printf("The heap is empty.\n");
        return;
    }

    match min = h->a[0];
    variant *v = hash_t[min.var]; // Ensure the index is valid and matches the variant

    if (v) {
        printf("Exact matches of any genetic variant were not found but the genome showcases a risk of:\n");
        printf("The variant of '%s'.\nIt was found to be located at %d position in the input genome.\nThe variant is of the type %s and attacks on chromosome number %d.\nNOTE: there is a high risk of this gene getting carried over to future generations!\n",
               v->gene, min.pos, v->significance, v->chromosome);
    } else {
       // printf("No valid variant found at the given index.\n");
       return;
    }
}


