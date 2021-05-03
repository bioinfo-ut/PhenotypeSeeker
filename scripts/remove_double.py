import itertools

kmer = "afdsfds"
kmer_coef = "43423"
samples_with_kmer = [4332, 434324]

a = "%s\t%s\t%s\t| %s\n" % (
                kmer, kmer_coef,
                2, " ".join(samples_with_kmer)
                )

print(a)