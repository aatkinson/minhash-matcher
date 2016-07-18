# MinHash Record Matcher

##### [For Sortable](http://sortable.com/challenge/)

### By Adam Atkinson

### To run: 
1. Have Python 2.6+ installed*.  
2. Clone this repo. 
3. Navigate to the root folder of this repo.
4. Run "go.sh".
5. Results are in "results.txt".

\* This implementation uses modules from the Python standard library only (i.e. everything included with the default Python installation).

### Notes

Minhashing is a scheme used to hash records quickly based on how similar they are. Essentially a record has a vector associated with it which encodes the words it contains. The minhash is the first non-zero element in a permutation of this vector. Minhashing is super cool because the probability that 2 records will cause a hash collision is equal to the Jaccard similarity of the records! You can use multiple permutations to build longer minhash signatures, and switch up the permutations you use to tweak the probability of having a hash collision when records have a certain similiarity.

These minhash signatures are computed for records and used as the basis of locality sensitive hashing (LSH) which matches the records together by similarity - which is made implicit by minhash.

Ok, that sounds wicked and all, but sounds overkill. Why'd I choose minhash?

1. The challenge involves approximate string matching on a lot of records. I felt that sophisticated fuzzy matching with edit distances would take too long in practice.

2. We need to avoid false positives more than false negatives. This means we need an algorithm that we can tweak to achieve this - and minhashing has this! With minhashing we have an S-curve defining the hash collision probability distribution which illustrates the false positive and negative prevalence. This S-curve can easily be defined with some parameters.

3. It mimics real life solutions to this problem. Here we pre-compute the hash table for the potential matches from a target data set and can match items quickly in constant time with hashing. Minhashing has offline extraction and online querying like many Big Data solutions. This technique also fits into the Map Reduce paradigm nicely.

4. I've been stoked to try out minhash since I first learned about it in my Big Data class at the University of Waterloo this past winter. (Check out the class here: lintool.github.io/bigdata-2016w/)

Minhash and LSH are awesome, Google them for some more details if you're curious!

Anyway, I'll talk about my implementation a bit here.

1. Read the products file to get the vocabulary of words we are looking for and give them integer IDs (vocab IDs - or 'vids') to make indexing easier. Also build a characteristic matrix to save on another read and tokenization.

2. Create the hash functions we'll use to permute the vectors representing the records. Also choose band and row parameters for the lookup table.

3. Read the product file again and compute the minhash signature for each line. Set up a the banded minhash lookup table for faster LSH and more precision.

4. For each entry in the listings file, compute its minhash signature and look it up in the banded minhash lookup table to see if there's any matches. 

5. Write the results to an output file.

#### Notes about threshold parameters:

You can change the characteristics of the S-curve with 2 command line args: sim-thresh and pr-@-sim-thresh
Basically, 'sim-thresh' is the Jaccard similarity and 'pr-@-sim-thresh' is the probability of a hash collision at 'sim-thresh.' Together these parameters characterize the S-curve.

By default these are:
- 0.99 for 'sim-thresh' because we really want the listing to have as many words in the desired product as possible (increased precision).
- 0.99 for 'pr-@-sim-thresh' because we really don't want to miss matches (increased recall).

If you want to loosen the precision, try values of 0.9, 0.8, ... , 0.5 for 'sim-thresh'.
   
#### Other notes:

- All done in standard Python 2.6+ to ease portability
- Chose Python because I'm good at it and it's fast to code things. I know it's slow.
- I normalize Unicode while parsing. This means that a letter with an accent just gets the accent removed.
- I don't have much experience with proper Unicode in Python. My implementation seems to work with Unicode though.
- I've tried to minimize the number of times I'm reading the input files. The nature of the required output file format necessitates more reading. 
