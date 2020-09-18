function DATASET = smallbone2011toyData()

DATASET.DESCRIPTION = 'smallbone2011 Toy test data.';

DATASET.FILE = 'smallbone2011_data';

DATASET.GROUPS = {
    'smallbone2011toy'
};

DATASET.FIELDS = {
    's1' 1 't' 's1_mean' 's1_std' 1 [] %Glc
    's2' 1 't' 's2_mean' 's2_std' 1 [] %G1P
    's3' 1 't' 's3_mean' 's3_std' 1 [] %G6P
    's4' 1 't' 's4_mean' 's4_std' 1 [] %Trehalose
    's5' 1 't' 's5_mean' 's5_std' 1 [] %T6P
    's6' 1 't' 's6_mean' 's6_std' 1 [] %UDP_glucose
};