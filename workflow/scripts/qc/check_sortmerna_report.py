#! /usr/local/env python3
import pandas as pd
import sys
import traceback



def main(exp_ids, logs, threshold, out_f):
    error_message = [""]
    merged_df = None
    for i, report in enumerate(logs):
        smr_df = pd.read_csv(report, header = 3, sep = '\t')
        smr_df.to_csv('/data/brinkvd/isabl/data_lake/analyses/27/20/52720/test2.csv')
        print(smr_df.columns)
        smr_df['exp_id'] = exp_ids[i]
        if smr_df['mean_perc'] > threshold:
            error_message += ["Error in sortmerna report for ",
                               exp_ids[i],
                               " value of sortmerna reads was: ",
                               smr_df['mean_perc'], '\n']
        if merged_df:
            merged_df = merged_df.append(smr_df, ignore_index=True)
        else:
            merged_df = smr_df
    merged_df.to_csv(out_f)
    with open(out_f, 'a') as o_f:
        o_f.write(''.join(error_message))


if __name__ == "__main__":
    
    with open(snakemake.log.e, "w") as ef, open(snakemake.log.o, "w") as of:
        sys.stderr = ef
        sys.stdout = of

        try:
            main(
                snakemake.params.experiment_ids,
                snakemake.input.smr_logs,
                snakemake.params.threshold,
                snakemake.output.sortmerna_report,
            )
        except Exception as e:
            traceback.print_exc(file=ef)