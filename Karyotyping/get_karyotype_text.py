import os
import pandas as pd


def sort_arms(segment_df):
    # Derive column 'arm' from column: 'band'
    # Transform based on the following examples:
    #    band        Output
    # 1: "p36.32" => "p"
    segment_df.insert(3, "arm", segment_df["band"].str[:1])
    # Derive column 'section' from column: 'band'
    # Transform based on the following examples:
    #    band        Output
    # 1: "p36.32" => "36.32"
    # segment_df.insert(3, "section", segment_df.apply(lambda row : row["band"][row["band"].find(".") - 2:], axis=1))
    segment_df.insert(3, "section", segment_df["band"].str[1:])
    # Sort by columns: 'seqnames' (ascending), 'arm' (ascending), 'section' (ascending)
    segment_df = segment_df.sort_values(['seqnames', 'arm', 'section'])
    # Change column type to float64 for column: 'section'
    # segment_df = segment_df.astype({'section': 'float64'})
    return segment_df


def parse_number(num_str):
    """Convert the series to a tuple (section, part) for easier comparison"""
    if '.' in num_str:
        section, part = num_str.split('.')
        return (int(section), int(part))
    else:
        return (int(num_str), None)  # Whole section, no part

def convert_to_string(series, arm):
    # Apply the parsing function
    parsed_series = series.apply(parse_number)

    # Sort the parsed series
    parsed_series = sorted(parsed_series)

    # Group consecutive numbers
    ranges = []
    start = parsed_series[0]
    for i in range(1, len(parsed_series)):
        # Check if they belong to the same section and are consecutive (with parts)
        if parsed_series[i][0] != parsed_series[i-1][0] or \
        (parsed_series[i][1] != parsed_series[i-1][1] + 1 if parsed_series[i][1] is not None else True):
            # If the current number is not consecutive, add the range and reset
            ranges.append((start, parsed_series[i-1]))
            start = parsed_series[i]
    # Add the last range
    ranges.append((start, parsed_series[-1]))

    # Format the ranges into the required string format
    def format_range(r):
        # eg. dup(5)(p14p15.3),dup(5)(p14p15.3)
        if r[0][1] is None and r[1][1] is None:  # Both are whole sections
            return f"({arm}{r[0][0]})"  # Only the section, no range
        elif r[0][0] == r[1][0] and r[0][1] == r[1][1]:  # Both parts are the same 
            return f"({arm}{r[0][0]}.{r[0][1]})"  # Only the section, no range
        elif r[0][1] is None:  # Only the start is a whole section
            return f"({arm}{r[0][0]}{arm}{r[1][0]}.{r[1][1]})"
        elif r[1][1] is None:  # Only the end is a whole section
            return f"({arm}{r[0][0]}.{r[0][1]}{arm}{r[1][0]})"
        else:  # Both are part sections
            return f"({arm}{r[0][0]}.{r[0][1]}{arm}{r[1][0]}.{r[1][1]})"

    # Join all ranges
    result = ''.join([format_range(r) for r in ranges])
    
    return result


def add_chrom_copies_to_string(result, chrom, copies):
    """
    if there is a double duplication of the short arm of chromosome 5 from bands p14 to p15.3, it would be written as:
    46,XX,dup(5)(p14p15.3),dup(5)(p14p15.3)
    """
    if copies == 1:
        result = f"del({chrom}){result},del({chrom}){result},del({chrom}){result}"
    if copies == 2:
        result = f"del({chrom}){result},del({chrom}){result}"
    if copies == 3:
        result = f"del({chrom}){result}"
    if copies == 5:
        result = f"dup({chrom}){result}"
    if copies == 6:
        result = f"dup({chrom}){result},dup({chrom}){result}"
    if copies == 7:
        result = f"dup({chrom}){result},dup({chrom}){result},dup({chrom}){result}"
    if copies == 8:
        result = f"dup({chrom}){result},dup({chrom}){result},dup({chrom}){result},dup({chrom}){result}"
    
    return result


name = "C3_Control"

fp = os.path.join(os.getcwd(), f'{name}_karyotype.perband.csv')
write_fp = os.path.join(os.getcwd(), f'{name}_karyotype.text.csv')
with open(fp, 'r') as fd:
    segment_df = pd.read_csv(fd, index_col=0)
    
segment_df_sorted = sort_arms(segment_df.copy())

# Grouped on columns: 'seqnames', 'Copies'
segment_df_grouped = segment_df_sorted.groupby(['seqnames', 'Copies', 'arm'])
fw = open(write_fp, "w")
fw.write(f"Complex Karyotype for {name},\n")
for df_group in segment_df_grouped:

    # print(df_group)
    chrom_num = df_group[0][0]
    copies = df_group[0][1]
    arm = df_group[0][2]
    # print(chrom_num, copies, arm)
    # print(type(df_group[1]['section']))
    # print(df_group[1]['section'])
    if copies == 4:
        continue
    
    # Generate the string
    result = convert_to_string(df_group[1]['section'], arm)

    result_w_chr = add_chrom_copies_to_string(result, chrom_num, copies)

    print(result_w_chr)
    fw.write(f"{result_w_chr},\n")
    
    # break
fw.close()
