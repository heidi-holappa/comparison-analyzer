
def create_dict_of_transcripts_and_reads(original_read_list: str, transcript_set: set):
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
                    dict_of_results[line.rstrip(
                        "\n").split("\t")[1]] = set()
                dict_of_results[line.rstrip("\n").split("\t")[1]].add(
                    line.rstrip("\n").split("\t")[0])
    return dict_of_results


def create_dict_of_reads_and_references(intron_site_dict: dict, dict_of_transcripts_and_reads: dict):
    reads_and_references = {}
    for intron_site_key, intron_site_values in intron_site_dict.items():
        if not intron_site_values['transcript_id'] in dict_of_transcripts_and_reads:
            continue
        for read in dict_of_transcripts_and_reads[intron_site_values['transcript_id']]:
            if read not in reads_and_references:
                reads_and_references[read] = set()
            reads_and_references[read].add(intron_site_key)

    return reads_and_references
