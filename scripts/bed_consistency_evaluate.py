import pandas as pd

def read_bed(filename):

    try:

        df = pd.read_csv(filename, sep='\t', header=None, usecols=[0, 1, 2], names=['chr', 'start', 'end'])
        return df
    except FileNotFoundError:
        print(f"Error: The file {filename} does not exist.")
        return None
    except pd.errors.EmptyDataError:
        print("Error: The file is empty.")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

tp = read_bed('/data/zjjiang/GRCH38_hificcl_pan_bed/80/tp.bed')
fp = read_bed('/data/zjjiang/GRCH38_hificcl_pan_bed/80/fp.bed')
fn = read_bed('/data/zjjiang/GRCH38_hificcl_pan_bed/80/fn.bed')


def filter_intervals(df, min_len, max_len):

    lengths = df['end'] - df['start']
    return df[(lengths >= min_len) & (lengths < max_len)]



ranges = [(500, 1000),(1000, 3000), (3000, 7000), (7000, float('inf'))]

for min_len, max_len in ranges:
    tp_filtered = filter_intervals(tp, min_len, max_len)
    fp_filtered = filter_intervals(fp, min_len, max_len)
    fn_filtered = filter_intervals(fn, min_len, max_len)

    # 计算数量
    TP = len(tp_filtered)
    FP = len(fp_filtered)
    FN = len(fn_filtered)


    precision = TP / (TP + FP) if (TP + FP) > 0 else 0


    recall = TP / (TP + FN) if (TP + FN) > 0 else 0

    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0


    print(f"Range {min_len}-{max_len}bp: Precision: {precision:.4f}, Recall: {recall:.4f}, F1 Score: {f1_score:.4f}")