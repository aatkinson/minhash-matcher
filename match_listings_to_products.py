#!/usr/bin/python
################################################################################
# match_listings_to_products.py
#
# A record linker based on minhash, implemented by Adam Atkinson as a solution
# to the Sortable coding challenge (http://sortable.com/challenge/).
#
# To run, open in a Terminal and: 
# `python match_listings_to_products.py products listings [sim-thresh] [pr-@-thresh]
#
# Produces an output file 'results.txt' with each line for one product and all
# of it's matched listings.
# 
# Two optional args are:
# - sim-thresh:  x-value of the S-curve which is the desired Jaccard similarity
# - pr-@-thresh: the probability of a hash collision at 'sim-thresh'
#
# See README for more details.
################################################################################

import codecs
import json
import linecache
import math    
import random
import re
import sys
import unicodedata

################################################################################
#   HELPER FUNCTIONS
################################################################################

def tokenize_raw_string (s, normalize_unicode=True, ignore_unicode=False):
    """
    Tokenizes a string by putting it to lowercase, normalizing unicode, and
    extracting alphanumeric substrings.

    Parameters:
        s: the raw string.
        normalize_unicode: whether or not to convert accented chars to unaccented/
        ignore_unicode: whether or not to skip unicode characters

    Returns:
        A list of ASCII string tokens.
    """
    s = s.lower()

    # This implementation works with ASCII only, choose what to do with unicode 
    if normalize_unicode:
        s = unicodedata.normalize('NFKD', s).encode('ascii','ignore')

    elif ignore_unicode:
        s = (unicode(s)).encode('ascii','ignore')

    return re.findall(r'[a-zA-Z0-9]+', s)

def get_prime_bigger_than_vocab (n):
    """
    Computes the smallest Mersenne prime that is larger than the input.

    Parameters:
        n: a non-negative integer.

    Returns:
        A (Mersenne) prime number.
    """
    x = math.ceil(math.log(n+1,2)) # Need prime >= |vocab|, n+1 prevents pow 2 case
    return int(math.pow(2,x)) - 1

def match_tokens_to_vids (tokens, vocabulary):
    """
    Matches words to word IDs (vids).

    Parameters:
        tokens: a list of ASCII strings.
        vocabulary: a dictionary of { "string":int} pairs representing the map.

    Returns:
        A list of ints (vids) corresponding to the tokens.
    """
    return [ vocabulary[vtok] for vtok in tokens if vtok in vocabulary]

def minhash_matched_vids (matched_vids, hfns):
    """
    Computes the minhash signature for some word IDs (vids).

    Parameters:
        matched_vids: a list of non-negative integers (vids).
        hfns: a list of hash functions which each take an integer.

    Returns:
        The minhash signature as a tuple indexed by the hash function used.
    """
    minhash_signature = [sys.maxint for h in hfns]

    for h, hfn in enumerate(hfns):
        for i in matched_vids:
            cur_hash = hfn(i)
            if cur_hash < minhash_signature[h]:
                minhash_signature[h] = cur_hash

    return tuple(minhash_signature)    

def minhash_tokenized_item (tokens, vocabulary, hfns):
    """
    Computes the minhash signature for a list of tokens.

    Parameters:
        tokens: a list of ASCII strings.
        vocabulary: a dictionary of { "string":int} pairs representing the map.
        hfns: a list of hash functions which each take an integer.

    Returns:
        The minhash signature as a tuple indexed by the hash function used.
    """
    return minhash_matched_vids(match_tokens_to_vids(tokens, vocabulary), hfns)

def find_matches_with_lsh (minhash, lookup_table, n_bands, n_rows):
    """
    Does locality sensitive hash to find similar items.

    Parameters:
        minhash: the minhash signature to lookup, as a tuple.
        lookup_table: a 2D dictionary keyed first by band, second by minhash part
        n_bands: number of bands in the lookup_table.
        n_rows: number of rows in the loopup_table (i.e. size of minhash part).

    Returns:
        A list of ints corresponding to the product IDs that are similar.
        (i.e. the lines in the products file that are similar)
    """
    similar_items = []
    for b in xrange(n_bands):
        minhash_section = tuple([ minhash[r] \
                                for r in xrange(b*n_rows, (b+1)*n_rows)])

        if minhash_section in lookup_table[b]:   
            new_similar_items = lookup_table[b][minhash_section]
            similar_items += new_similar_items
        else:
            continue

    # De-duplicate the matched items and return as a list
    return list(set(similar_items))

################################################################################
#   MAIN PROGRAM
################################################################################

if __name__ == "__main__":

    random.seed()
    
    if len(sys.argv) < 3:
        print "ERROR: not enough input arguments!"
        print "Usage: python match_listings_to_products.py products listings [0<#<1 similarity threshold] [0<#<1 probability of detection]"
        sys.exit(1)

    prod_file = sys.argv[1]
    listings_file = sys.argv[2]
    results_file = "results.txt"

################################################################################
#   Part 1: Get vocabulary list and inverted characteristic matrix
################################################################################

    n_products = 0
    vocabulary = {}
    vocab_size = 0
    product_names = []

    with codecs.open(prod_file, 'r', 'utf-8') as f:

        # First pass through products to get the size of the vocabulary
        # and to save the product names. 
        # We'll need the product names to generate the output, and better to
        # keep them in memory than the much larger listings file
        for i, line in enumerate(f):
            line_object = json.loads(line)
            product_names.append(line_object['product_name'])
            n_products += 1

        matched_vids_map = [ None for p in xrange(n_products)]

        # Next parse lines for characteristic matrix entries
        # We'll do this here to save an extra read & parse of the products file
        f.seek(0)
        for i, line in enumerate(f):
            line_object = json.loads(line)

            vocab_tokens = tokenize_raw_string(line_object['product_name'])
            vocab_tokens += tokenize_raw_string(line_object['manufacturer'])
            vocab_tokens += tokenize_raw_string(line_object['model'])

            vids = []
            for vtok in vocab_tokens:

                # If key is new, update the maps
                if vtok not in vocabulary:
                    vocabulary[vtok] = vocab_size
                    vocab_size += 1

                vids.append(vocabulary[vtok]) # dupes possible

            matched_vids_map[i] = vids

################################################################################
#   Part 2: Set up the hash functions
################################################################################

    # Get the desired Jaccard similarity so we can set choose our N/bands and
    # K/rows to set the S curve threshold
    similarity_threshold = 0.975
    pr_at_similarity_threshold = 0.99
    if len(sys.argv) >= 5:
        try:
            similarity_threshold = float(sys.argv[3])
            pr_at_similarity_threshold = float(sys.argv[4])
        except ValueError as e:
            print "ERROR: invalid thresholds!"
            sys.exit(1)

    if not ((similarity_threshold > 0 and similarity_threshold < 1) and \
            (pr_at_similarity_threshold > 0 and pr_at_similarity_threshold < 1)):
        print "ERROR: invalid thresholds!"
        sys.exit(1)

    # Calculate b & r from p = 1-(1-s^r)^b
    # ***NOTE: if you want to enter # of bands and rows manually, do it here***
    r_rows = 10
    b_bands = int(math.ceil(math.log(1-pr_at_similarity_threshold) \
                / math.log(1-math.pow(similarity_threshold, r_rows))))
    n_hash_fns = r_rows * b_bands

    # Get the prime number for our hash functions
    prime = get_prime_bigger_than_vocab(vocab_size)

    print "\nMinhashing Parameters:\n"
    print "Total number of hash functions for minhash signatures: {0}".format(n_hash_fns)
    print "Jaccard similarity threshold: {0:.2f}%".format(100*similarity_threshold)
    print "Probability that a document with Jaccard similarity is detected: {0:.2f}%"\
                                        .format(100*pr_at_similarity_threshold)
    print "Number of bands in lookup table (b): {0}".format(b_bands)
    print "Number of rows in each band of the lookup table (r): {0}".format(r_rows)

    # Hash function generating function
    def uniform_hash (a,b,p,m):
        return lambda x: (((a*x) + b) % p) % m

    hash_functions = [  uniform_hash(random.randint(1,sys.maxint), \
                        random.randint(0,sys.maxint),prime, vocab_size) \
                        for k in xrange(n_hash_fns)]

################################################################################
#   Part 3: Minhash all of the products
################################################################################

    print "\nBuilding minhash signatures and lookup table for products...\n"

    minhash_signature_mat = [ [sys.maxint]*n_products for h in hash_functions ]
    banded_minhash_lookup_map = [ {} for b in xrange(b_bands) ]

    for p in xrange(n_products):

        # Compute the minhash signature and insert into the signature matrix
        minhash_sig = minhash_matched_vids(matched_vids_map[p], hash_functions)

        for h in xrange(len(hash_functions)):
            minhash_signature_mat[h][p] = minhash_sig[h]

        # Speedup product lookup by hashing the minhashes (in a Python dict) 
        for b in xrange(b_bands):
            for r in xrange(r_rows):

                signature_section = tuple([ minhash_signature_mat[r][p] \
                                    for r in xrange(b*r_rows, (b+1)*r_rows)])

                if signature_section not in banded_minhash_lookup_map[b]:
                    banded_minhash_lookup_map[b][signature_section] = []

                banded_minhash_lookup_map[b][signature_section].append(p)

        if p % 100 == 0:
            print "> {0}/{1} products processed".format(p, n_products)

################################################################################
#   Part 4: Match the listings with locality sensitive hashing
################################################################################

    print "\nMatching listings against products...\n"

    # Build the inverted product:matched_listings map, 
    # compared to listings:matched_products
    product_to_matched_listing_indices = [ [] for i in xrange(n_products)]

    listings_matched = 0
    listings_processed = 0

    # Read each listing and minhash it against the vocabulary from the products
    with codecs.open(listings_file, 'r', 'utf-8') as f:

        for i, line in enumerate(f):
            line_object = json.loads(line)
            vocab_tokens = tokenize_raw_string(line_object['title'])
            vocab_tokens += tokenize_raw_string(line_object['manufacturer'])

            item_minhash = minhash_tokenized_item(  vocab_tokens, \
                                                    vocabulary, \
                                                    hash_functions )
            similar_items = find_matches_with_lsh(  item_minhash, \
                                                    banded_minhash_lookup_map, \
                                                    b_bands, \
                                                    r_rows )

            # Break multiple product matches randomly
            if similar_items:
                matched_product = random.choice(list(similar_items))
                product_to_matched_listing_indices[matched_product].append(i)
                listings_matched += 1

            listings_processed += 1

################################################################################
#   Part 5: Generate output file
################################################################################

    skip_unmatched_products = False
    
    linecache.clearcache()

    # Get the listings by random access to the file via linecache module
    with codecs.open(results_file, 'w', 'utf-8') as outfile:

        json_encoder = json.JSONEncoder(ensure_ascii=False)

        for p in xrange(n_products):

            if skip_unmatched_products and not product_to_matched_listing_indices[p]:
                continue

            linecache.checkcache(listings_file)
            result = { 'product_name':product_names[p], 'listings':[] }
            
            for i in product_to_matched_listing_indices[p]:
                listing_json = json.loads(linecache.getline(listings_file, i+1))
                result['listings'].append(listing_json)

            result_json_str = json_encoder.encode(result)
            outfile.write(result_json_str+"\n")

    print "Total percentage of listings matched = {0:.2f}%\n".format(100*float(listings_matched)/max(1,listings_processed))
