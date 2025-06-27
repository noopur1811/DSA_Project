🧬 DNA Variant Detection using Trie and Aho-Corasick Algorithm

This project implements a **genome variant detection system** using **advanced Data Structures and Algorithms** in C, primarily using **Trie**, **Aho-Corasick pattern matching**, **hash tables**, and **heaps** for approximate matching. It can detect **exact or closest variants** in DNA sequences, and print detailed reports for individual users based on their genomic sequence.


````markdown
## 📁 Project Structure

```bash
.
├── main.c                  # Driver code
├── trie.c / trie.h         # Trie implementation for variant sequences
├── heap.c / heap.h         # Min-heap for storing best approximate matches
├── hash.c                  # Hashing for storing reference genome metadata
├── filled_reference_genome.csv # List of labeled known gene variants
├── trial.csv              # User test sequences with name, age, gender
├── a.exe                  # Compiled executable (Windows)
````

---

## 🚀 Features

* 📌 **Exact Match Detection** using Aho-Corasick on DNA sequences (A, C, G, T)
* 💡 **Approximate Matching** using edit distance for near matches
* 🧠 **Heap** stores top match with minimum edit distance
* 📋 **Hash Table** stores metadata like frequency, label, significance of each variant
* 🧾 **Personalized Report Generation** for each user
* 📤 **Command-Line Interface** output with rich formatting

---

## 🧪 Sample Output

```
---------------------------------------------------------------------------------------------
USER: 8
 Name: Michael Ray
 Age: 45
 Gender: Male

Report:
Identified Gene is of the variant name: Down Syndrome
Sequence of chromosomes in the variant was found to be: AATGTTAAAGTTCTAGTGGA
The gene is of Pathogenic type
Frequency of the identified variant is 1.20
Position at which the genome is affected is 79
---------------------------------------------------------------------------------------------

USER: 9
 Name: Laura King
 Age: 31
 Gender: Female

Report:
Exact matches of any genetic variant were not found, but the genome showcases a risk of:
The variant of 'BRCA2'.
It was found to be located at position 41 in the input genome.
The variant is of the type Pathogenic and attacks on chromosome number 9.
NOTE: there is a high risk of this gene getting carried over to future generations!
---------------------------------------------------------------------------------------------
```

---


---

## ⚙️ How It Works

1. 🧱 **Trie Construction**  
   All known gene variant sequences are inserted into a trie.

2. 🔁 **Aho-Corasick Failure Links**  
   Built for fast multi-pattern detection with backtracking.

3. 🔬 **User Sequence Scanning**  
   Each test sequence is matched against the trie for exact or near matches.

4. 🧮 **Match Handling**
   - ✅ **Exact Match** → Report immediately
   - ❌ **No Match** → Approximate match via edit distance + heap

5. 🧾 **Report Generation**  
   Based on variant metadata from the hash table.

---

## 📊 Input File Formats

### `filled_reference_genome.csv`


## 📊 Input Files

* `filled_reference_genome.csv`
  Format: `Key,Label,Sequence,Significance,Frequency`

* `trial.csv`
  Format: `Name,Age,Gender,Test DNA Sequence`

---

## 🛠️ Build & Run

```bash
# Compile
gcc main.c trie.c heap.c hash.c -o dna_variant

# Run
./dna_variant
```

> Or double-click `a.exe` on Windows to run directly.

---

## 👨‍💻 Contributors
* **Ankita Karhade** — 
* **Noopur Karkare** — [GitHub](https://github.com/noopur1811)

---

📚 Topics Covered
Aho-Corasick Pattern Matching

Edit Distance (Levenshtein)

Trie Construction and Search

Hash Tables and Metadata Lookup

Min-Heap for Closest Match

File Handling and CSV Parsing in C

CLI-based Report Generation


---

## 📄 License

This project is for educational purposes and is licensed under the MIT License.

```
