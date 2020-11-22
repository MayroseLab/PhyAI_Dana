import numpy as np
import pandas as pd
#import pandas_profiling as pp
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import itertools

from sklearn.model_selection import KFold, RandomizedSearchCV
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn import metrics

from scipy.sparse import csr_matrix

pd.set_option('display.max_columns', None)
sns.set_style("white")
palette = itertools.cycle(sns.color_palette('colorblind'))

import warnings
warnings.filterwarnings('ignore')

DIRPATH = r'C:\Users\ItayMNB3\Desktop\Assignment_Microsoft\data\\'
VARS_DICT = {'bdate': 'DateOfBirth', 'subj': 'SubjectId', 'name': 'Name', 'parent': 'ParentId', 'level': 'Level',
             'ques':'QuestionId', 'user':'UserId', 'ans':'AnswerId', 'adate':'DateAnswered'}


def print_val_type(dataframe, colname, nrows=100):
    ''' get the type of list/list-like(str) values, as dtype will return obj (str) in any case '''
    for l in dataframe.head(nrows)[colname]:
        print("Type of col {} is:".format(colname), type(l))

    # return the last occurance
    return type(l)


def return_sorted_lst(l, func):
    '''
    sort the elements of an IDsList by the tree structure (that shoud be provided as a dictionary to the input func)-
    *for a single-path list (namely exactly 4 elements: return a copy of the sorted list
    *for a path the stops in an internat node (namely subjectxx-other): return a copy of a list without the halting node (redundant information)
    *for a not-single path (>4 elements): return only the root node (just because! time-limitated task & <10% of the samples..)
    '''
    l.sort(key=func)
    treedeep = len(l)
    new_l = l[:-1] if treedeep < 4 else l[:1] if treedeep > 4 else l[:]

    return new_l


def remove_from_lst(lst_remove_from, lst_remove_these):
    ''' remove specified list elements from a given list'''
    for elem in lst_remove_these:
        if elem in lst_remove_from:
            lst_remove_from.remove(elem)

    return lst_remove_from


df_train = pd.read_csv(DIRPATH + r'train_data\train_task_1_2.csv')
df_meta_student = pd.read_csv(DIRPATH + r'metadata\student_metadata_task_1_2.csv')
df_meta_question = pd.read_csv(DIRPATH + r'metadata\question_metadata_task_1_2.csv')
df_meta_subject = pd.read_csv(DIRPATH + r'metadata\subject_metadata.csv')
df_meta_answer = pd.read_csv(DIRPATH + r'metadata\answer_metadata_task_1_2.csv')



# print some information
print(df_meta_question)
type_sample1 = print_val_type(df_meta_question, VARS_DICT['subj'], nrows=1)

# convert the value type to lists
if not type_sample1 == list:
    print("number of unique paths:",len(df_meta_question[VARS_DICT['subj']].unique()))
    df_meta_question[VARS_DICT['subj']] = df_meta_question[VARS_DICT['subj']].apply(eval)
    print("\nAfter converting-")
    print_val_type(df_meta_question, VARS_DICT['subj'], nrows=1)

    # validate that the proportion of the non-single-path lists is low, so I can make an "easy" (=quick) decision of how to handle it
    df_question2subj_path4 = df_meta_question[df_meta_question[VARS_DICT['subj']].map(len) > 4]
    print("\nThe percentages of samples that contain more than 4 features, namely more than one path, is:")
    print(round((100*len(df_question2subj_path4)/len(df_meta_question)),2), "%")



### I can easily tell that the lists are NOT "sorted" by its level.
### so I must go over each list and sort it (see ducumentation of 'return_sorted_lst' function above)

d_subj2level = pd.Series(df_meta_subject[VARS_DICT['level']].values, index=df_meta_subject[VARS_DICT['subj']]).to_dict()
def myFunc(e, d=d_subj2level):
    return d[e]

df_meta_question[VARS_DICT['subj']] = df_meta_question.apply(lambda x: return_sorted_lst(x[VARS_DICT['subj']], myFunc), axis=1)


### expand to cols and index subj_ids in the proper SubjID_deep_x col

res_df_questions = df_meta_question[VARS_DICT['subj']].apply(pd.Series).astype('category')
res_df_questions.insert(0, VARS_DICT['ques'], df_meta_question[VARS_DICT['ques']])
res_df_questions = res_df_questions.rename({i: 'SubjID_deep{}'.format(i) for i in range(4)}, axis='columns')




if VARS_DICT['bdate'] in df_meta_student.columns:            # namely if it is the first time I run this section
    df_meta_student['temp'] = pd.to_datetime(df_meta_student[VARS_DICT['bdate']], dayfirst=True, errors='coerce')
    df_meta_student['Byear'] = df_meta_student['temp'].dt.year

    #filtering illegal values (according to the pdf: questions provided between 2018 to 2020. students roughly between 7 and 18 years old
    df_meta_student.loc[df_meta_student['Byear'] < 2000] = None
    df_meta_student.loc[df_meta_student['Byear'] > 2012] = None

    df_meta_student['Bmonth'] = df_meta_student['temp'].dt.month.astype('category')

    df_meta_student = df_meta_student.drop(['temp', VARS_DICT['bdate']], axis=1)



if VARS_DICT['adate'] in df_meta_answer.columns:            # namely if it is the first time I run this section
    df_meta_answer['temp'] = pd.to_datetime(df_meta_answer[VARS_DICT['adate']])
    df_meta_answer['Ayear'] = df_meta_answer['temp'].dt.year.astype('category')
    df_meta_answer['Amonth'] = df_meta_answer['temp'].dt.month.astype('category')
    df_meta_answer['Day'] = df_meta_answer['temp'].dt.dayofweek.astype('category')

    ## making up (because of time limitation) a conversion dictionary to assign the 24hours to 4 time_period groups.
    ## I am sure there are solutions that make much more sense than the following one:
    DayTime_dict = {i:k for i,k in enumerate(np.concatenate((np.zeros(6), np.ones(6), np.full_like(np.arange(6, dtype=int), 2), np.full_like(np.arange(6, dtype=int), 3))))}
    df_meta_answer['DayTime'] = df_meta_answer['temp'].dt.hour
    ## the following line takes a while- I should replace it with more efficient implementation (map? np.where? read) if I have the time when I finish
    df_meta_answer['DayTime'] = df_meta_answer.apply(lambda x: DayTime_dict[x['DayTime']], axis=1).astype('category')

    df_meta_answer = df_meta_answer.drop(['temp', VARS_DICT['adate']], axis=1)



train_data_len = len(df_train)
complete_df = df_train.merge(df_meta_student, on=VARS_DICT['user'], how='left').merge(df_meta_answer, on=VARS_DICT['ans'], how='left').merge(res_df_questions,on=VARS_DICT['ques'], how='left')

# I can add a feature: age_when_answering. return if I have more time


assert(train_data_len == len(complete_df))

# drop AnswerId because is unique and CorrectAnswer because..
features_to_drop = ['AnswerId', 'CorrectAnswer']
if features_to_drop[0] in complete_df.columns:
    complete_df = complete_df.drop(features_to_drop, axis=1)


# Actually all features but 'Confidence' are categorical. So convert thier type (in any case I will one-hot encode them):
cols_lst = list(complete_df)
for e in cols_lst:
    if e == 'Confidence':
        continue
    complete_df[e] = complete_df[e].astype('category')



special_cols = ['QuestionId', 'UserId', 'IsCorrect']
cols_high_card = ['GroupId','QuizId','SchemeOfWorkId','SubjID_deep2','SubjID_deep3']
label = special_cols[-1]
pivot_cols = special_cols[:-1]

# deal with features with extreme cardinality
#for col in cols_high_card:
#    print(col, len(complete_df[col].unique()))
########################### removing features with high card - TEMP! ###########################
## After implementing a base model: expand each of these to few most common values, and 'others'.
if cols_high_card[0] in complete_df.columns:
    complete_df = complete_df.drop(cols_high_card, axis=1)
################################################################################################

features_withCat = remove_from_lst(list(complete_df), special_cols)
complete_df_encoded = pd.get_dummies(complete_df, columns=features_withCat, prefix=features_withCat, sparse=True).fillna(0)
features = remove_from_lst(list(complete_df_encoded), special_cols)




## check memory usage after action
#complete_df_encoded.info(memory_usage='deep')

########## cause memory error - int overflow #########
df_transformed_X = complete_df_encoded.head(1000).pivot(
    index=pivot_cols[1],
    columns=pivot_cols[0],
    values=features)#.fillna(0)

df_transformed_Y = complete_df_encoded.head(1000).pivot(
    index=pivot_cols[1],
    columns=pivot_cols[0],
    values=[label])#.fillna(0)
####################################################


########## additional attempts #########
#df_transformed_x = complete_df_encoded.head(100).groupby([pivot_cols[1], pivot_cols[0]])[features].max().unstack()
#df_transformed_Y = complete_df_encoded.head(100).groupby([pivot_cols[1], pivot_cols[0]])[label].max().unstack()

#chunk_size = 1000
#chunks = [x for x in range(0, complete_df_encoded.shape[0], chunk_size)]
#df_transformed_X = pd.concat([complete_df_encoded.iloc[ chunks[i]:chunks[i + 1] - 1 ].pivot(
#    index=pivot_cols[1], columns=pivot_cols[0], values=features) for i in range(0, len(chunks) - 1)])
#df_transformed_Y = pd.concat([complete_df_encoded.iloc[ chunks[i]:chunks[i + 1] - 1 ].pivot(
#    index=pivot_cols[1], columns=pivot_cols[0], values=[label]) for i in range(0, len(chunks) - 1)])
########################################


# convert using the more efficient method of scipy sparse matrix
X_mat = csr_matrix(df_transformed_X.values)
Y_mat = csr_matrix(df_transformed_Y.values)


# Split dataset into training, test, and validation sets (80%, 10%, 10%)
X_train, X_test, y_train, y_test = train_test_split(df_transformed_X, df_transformed_Y, test_size=0.2, random_state=1)
X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.125, random_state=1) # 0.125 x 0.8 = 0.1