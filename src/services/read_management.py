def create_dict_of_transcripts_and_reads(transcript_set: set, original_read_list: str):
    """
    Creates a dictionary of read ids as keys and the transcripts they are mapped to.
    The relationship can be one-to-many, i.e. one read can be mapped to multiple transcripts.

    Args:
        output_filename_dict (dict): a dict with transcript ids as keys, output filenames as values
        original_read_list (str): a filename of the original read list

    Returns:
        dict: outputs a dictionary of read ids as keys and the transcripts they are mapped to
    """
    dict_of_results = {}
    with open(original_read_list, encoding="UTF-8") as file:
        for line in file:
            if line.rstrip("\n").split("\t")[1] in transcript_set:
                if line.rstrip("\n").split("\t")[1] not in dict_of_results:
                    dict_of_results[line.rstrip("\n").split("\t")[1]] = set()
                dict_of_results[line.rstrip("\n").split("\t")[1]].add(
                    line.rstrip("\n").split("\t")[0])
    return dict_of_results
