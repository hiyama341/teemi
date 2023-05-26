
# Importing the module we are  testing
from teemi.test.data_wrangling import *

sequencing_plates = pd.read_csv('../teemi/tests/files_for_testing/Plate2Seq.csv',index_col=False)
df = [sequencing_plates]

#Plate2SeqFunctions
sliced = slicing_and_naming_seq_plates(df)
list_of_dfs = plat_seq_data_wrangler(sliced)
filtered_plates = plate_AvgQual(list_of_dfs)
split_df = split_df_names(filtered_plates)
all_data_frames = concatenating_list_of_dfs(split_df)

def test_slicing_and_naming_seq_plates():
    
    assert len(sliced[0]) == 81


def test_plat_seq_data_wrangler():

    assert len(list_of_dfs[0]) == 81
    assert type(list_of_dfs[0].iloc[3]['AvgQual']) == type(np.float64(0))


def test_plate_AvgQual():

    # do we actually filter on our parameters?
    true_false = (filtered_plates[0].iloc[:]['AvgQual']>=50).any()
    true_false1 = (filtered_plates[0].iloc[:]['used']>=25).any()

    assert len(filtered_plates[0]) == 72
    assert true_false == True
    assert true_false1 == True


def test_split_df_names():
    assert  split_df[0].columns[7] == 'plate'
    assert  split_df[0].columns[8] == 'well'


def test_concatenating_list_of_dfs():
    assert len(all_data_frames) == 72
    assert len(all_data_frames.columns) == 9