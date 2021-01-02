import pandas as pd
pd.set_option('display.max_columns', 40)



#df = pd.read_csv(r"C:\Users\ItayM3\Dropbox\Anna.csv")
df = pd.read_csv("/groups/itay_mayrose/annarice/crispr-il/results/human/processed_for_accessibility/results/all_gRNA_human1000.bed", sep='\t', header=None)
desired_len = len(df[1].unique())

names_col6 = df[6].unique()
for new_col in names_col6:
    df[new_col] = df[df[6] == new_col].apply(lambda x:x[7], axis=1)

df = df.groupby(1).apply(lambda x: x.ffill().bfill()).drop_duplicates(subset=[1])
print(df.head(20))
assert desired_len == len(df)











