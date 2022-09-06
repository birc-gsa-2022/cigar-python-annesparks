# Edits and CIGARs

If we have strings `p` and `p'`, a string of "edits" can describe how to go from one to the other. For example, with `p = aa` and `p' = aba`, we could say that we can get from `p` to `p'` by keeping the first `a` in `p`, then inserting a `b`, and then keeping the second `a` (which is now the third). We would express this as a "match/mismatch" operation `M`, then an insert operation `I`, and then another match `M`:

```
     MIM
p  = a-a
p' = aba
```

If we only have `p` and the operations `MIM`, we cannot determine what `p'` is. Obviously, we don't know which letter we inserted in the `I` operation, but the `M` can both match and mismatch, so from `p` via `MIM` we could equally well have ended up with `bba`, `aca` or `aaa` and many more.

However, if we know that `p'` is a sub-string of a string `x`, and we know that it starts at index `i`, then the edit operations will tell us how long a prefix of `x[i:]` the transformation of `p` would match.

In this exercise, we will write some utility code that can help us work with approximative matches of `p` against `x`.

## Pairwise alignments

We will start with having strings `p` and `q` where `q` is a transformation of `p`. We can write this transformation in two different ways: as a pairwise aligment similar to above:

```
p = aacgt-c
q = a-attac
```

where we have two strings `p = aacgt-c` and `q = a-attac` of equal length, and where gap-characters `-` tells us where there are insertions or deletions, or we can represent the transformation as three strings

```
p = aacgtc
q = aattac
e = MDMMMIM
```

where the `e` string describes the edits to go from `p` to `q`.

Your first task is to write two functions that translates between the two representations.

```
get_edits("aacgt-c", "a-attac") -> ("aacgtc", "aattac", "MDMMMIM")
align("aacgtc", "aattac", "MDMMMIM") -> ("aacgt-c", "a-attac")
```

A small state-machine will get the job done for both functions. You look at the driver(s) of the machine, in the first function it is the front of the two strings and in the second it is the edits string, and then, depending on what that state is, you make a transformation and update the emitted strings.

For the first function, the transformations could look like this:

    1. If the two strings are empty, you are done.
    2. If the two strings both have non-gaps at the front, add the letters there to the
       corresponding output and put an `M` in the edits string.
    3. If the first string has a gap, then the other doesn't (that is an invariant we will
       insist on), so you add the second's string's letter to its corresponding output, and you
       add an `I` to the edits.
    4. If the second string has a gap, then add the first string's letter to its output and add
       a 'D' to the edits.

For the second function, the transformations could look like this:

    1. If there are no more edits, you are done.
    2. If you have an `M` edit, you emit the front letters from the two strings to the corresponding
       output.
    3. If you have an `I` edit, you emit a `-` to the first string's output and the second string's letter
       to the other output.
    4. If you have a `D` edit, you copy the first string's letter to its output and you emit a `-`
       for the second string.


You can also think of the functions as functional transformations. In the first function you can think of the state at the front of the two strings and apply the rule:

```haskell
get_edits (p, q, res_p, res_q, edits) = 
    case (p, q) of
        -- When we are through p and q (and we will at the same time), then
        -- we have our result (although we have built it in reverse here)
        ([], []) -> (rev res_p, rev res_q, rev edits),
        -- We can't have deletions in both strings at the same time (that is an
        -- invariant we will insist upon), so these are the remaining operations
        ( a  :: p', '-' :: q') -> get_edits(p', q', a :: res_p,      res_q, 'D' :: edits)
        ('-' :: p',  b  :: q') -> get_edits(p', q',      res_p, b :: res_q, 'I' :: edits),
        ( a  :: p',  b  :: q') -> get_edits(p', q', a :: res_p, b :: res_q, 'M' :: edits)
```

This is in "pseudo-Haskell" and you don't necessarily want to implement it this way--it could be rather inefficient in some languages--but the rules will be the same: you are done when you are through the strings, otherwise you consider the first letter in the two, handle a gap in the second as a deletion, a gap in the first as an insertion, and if you have non-gaps in both strings you have a match. Then you continue with the rest of the strings.

In the second function, the driver of the state-machine is the edits sequence, and the transformations can look like this:

```haskell
align (p, q, e, align_p, aligh_q) =
    case (e, p, q) of
        -- We are naturally done when there is nothing left and that
        -- should be e is empty. Something when wrong if either p or q
        -- is non-empty at this point
        (_, _, []) -> (rev align_p, rev align_q),
        -- Otherwise, the first letter in e tells us what to do
        ('M' :: e', a :: p', b :: q') -> align(p', q', e',  a  :: align_p,  b  :: align_q),
        ('I' :: e',      p', b :: q') -> align(p', q', e', '-' :: align_p,  b  :: align_q),
        ('D' :: e', a :: p',      q') -> align(p', q', e',  a  :: align_p, '-' :: align_q)
```


## Local alignments

In a read mapper, you don't have a complete (also known as "global") alignment, as the ones we handled above. Instead, we have an alignment that includes the entirety of `p` and a sub-string of `x`, `x[i:j]`. We typically do not know what the substring is, only where it starts, `x[i:?]`, but the string of edits contains enough information to work out the full string. If you ran the `get_edits` function from above on `p` and the suffix of `x` starting at `i`, `x[i:]`, and you didn't abort or threw an exception if `x` wasn't empty when you made it through the edits (as I didn't in the Haskell code), then whatever you had produced for the second string would tell you which substring `p` would be transformed into. Just remove the gaps from it, and you would have `x[i:j]`.

We won't bother with removing the gaps from the result, but for testing and validation purposes, it is useful to be able to extract the local alignment from the string `p`, the position `i`, and the string of edits, `e`.

Write a function `local_align(p, x, i, e) -> (p',x[i:j])` that gives you the local alignment of `p` against `x` starting at position `i` and following the edits in `e`.

I have found such a function quite useful for debugging purposes in the past, and you might as well.


## Edit distance

When searching for approximative matches, we often put restrictions on how far you can modify `p` to get `p'` and still call it a match--we certainly will in this class. If we want to automtically test that our algorithm is doing the right thing, it can be useful to validate that all the matches it finds are indeed within the required distance from the pattern. Luckily, it is straightforward to get the distance between `p` and `p'` if we have a local alignment. For each column, i.e. letter on top of another letter, count what that column contributes to the distance. If you have an insertion or a deletion, i.e., there is a `-` in one of the rows, you have used one edit. If you have two letters, then it isn't a true *edit* if they are the same, so that counts as zero, but if the two letters are different, then it is a substitution, and that counts.

The list of edits itself doesn't tell you how far you need to modify `p`, because you cannot distinguish between the free `M` operations that are matches and the mismatch `M`s that needs to be counted, but if you use your `local_align(p, x, i, e)` function from above, and then run through the two resulting strings, it is easy enough. Write a function, `edit_dist(p, x, i, e)`, that counts the number of edits.

With this function in hand, you can easily check that all the matches your algorithm finds, much later in the class, are truly matches.


## Run-length encoding

Read mappers do not use strings of edits the way we have done. With long reads and few edits, these are waseful in memory, and the output files are pretty large as it is. Instead, the edit operations are run-lenght encoded as CIGAR strings.

You will also have to use CIGAR strings in the final project, but edits are much easier to work with, so we want some functions to translate between the two representations. Write two functions, `edits_to_cigar` and `cigar_to_edits` that translate between the two representations:

```
edits_to_cigar("MMIMMDDM") -> "2M1I2M2D1M"
cigar_to_edits("2M1I2M2D1M") -> "MMIMMDDM"
```


