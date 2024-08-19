# HE_ReLU
For CSNOVA, 2024 summer, a program using the HE lib troy-nova to encrypt ReLU function.

## Next step
    * I have basically gotten the hang of how to utilize troy for CKKS, next job is to link Polynomial with CKKS.
    * The current way of poly-approximation is using Taylor's equation, but it does not have a very good effect. I may consider using Remez algorithm instead, which I haven't learnt how to use it yet.

## Some questions to think about
    * How to form the modular chain? According to the experiment we've conducted, I found that using 60 followed by a lot of 40s is effective for the fixed polynomial, but why? 
    * The polynomial fixed is ok when its deg reach 6, but what about the automated one? 