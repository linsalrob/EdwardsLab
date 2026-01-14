"""

"""

import os
import sys
import re
import argparse
import logging
__author__ = 'Rob Edwards'

from transformers import T5Tokenizer, AutoModelForSeq2SeqLM
import torch

logging.basicConfig(level=logging.INFO)


def run(args):
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    # Load the tokenizer
    tokenizer = T5Tokenizer.from_pretrained(args.model, do_lower_case=False)

    # Load the model
    model = AutoModelForSeq2SeqLM.from_pretrained(args.model).to(device)
    logging.info(f"Loaded model {args.model} on device {device}")
    logging.info(f"Model config: {model.config}")
    logging.info(f"Model is encoder/decoder: {model.config.is_encoder_decoder}")


    # only GPUs support half-precision currently; if you want to run on 
    # CPU use full-precision (not recommended, much slower)
    model.full() if device=='cpu' else model.half()

    # prepare your protein sequences/structures as a list.
    # Amino acid sequences are expected to be upper-case ("PRTEINO" below)
    # while 3Di-sequences need to be lower-case.
    sequence_examples = ["PRTEINO", "SEQWENCE"]
    min_len = min([ len(s) for s in sequence_examples])
    max_len = max([ len(s) for s in sequence_examples])

    # replace all rare/ambiguous amino acids by X (3Di sequences does not 
    # have those) and introduce white-space between all sequences (AAs and 3Di)
    sequence_examples = [" ".join(list(re.sub(r"[UZOB]", "X", sequence))) 
                         for sequence in sequence_examples]

    # add pre-fixes accordingly. 
    # For the translation from AAs to 3Di, you need to prepend "<AA2fold>"
    sequence_examples = [ "<AA2fold>" + " " + s for s in sequence_examples]

    # tokenize sequences and pad up to the longest sequence in the batch
    ids = tokenizer.batch_encode_plus(sequence_examples,
                                      add_special_tokens=True,
                                      padding="longest",
                                      return_tensors='pt').to(device)

    # Generation configuration for "folding" (AA-->3Di)
    gen_kwargs_aa2fold = {
                      "do_sample": True,
                      "num_beams": 3, 
                      "top_p" : 0.95, 
                      "temperature" : 1.2, 
                      "top_k" : 6,
                      "repetition_penalty" : 1.2,
    }

    # translate from AA to 3Di (AA-->3Di)
    with torch.no_grad():
        translations = model.generate( 
            ids.input_ids, 
            attention_mask=ids.attention_mask, 
            max_length=max_len, # max length of generated text
            min_length=min_len, # minimum length of the generated text
            early_stopping=True, # stop early if end-of-text token is generated
            num_return_sequences=1, # return only a single sequence
            **gen_kwargs_aa2fold
        )
        # Decode and remove white-spaces between tokens
    decoded_translations = tokenizer.batch_decode(translations, skip_special_tokens=True)
    # predicted 3Di strings 
    structure_sequences = [ "".join(ts.split(" ")) for ts in decoded_translations ]
    for i,s in enumerate(sequence_examples):
        print(f"Sequence: {s.replace(' ','')} 3Di: {structure_sequences[i].replace(' ','')}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-m', '--model', help='model', default="Rostlab/ProstT5")
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    run(args)


