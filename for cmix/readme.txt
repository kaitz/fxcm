This is model for cmix-hp for replacing paq8hp.

Notes:
In predictor.c interval6 map should point to wrt_t and interval4 wrt_w in fxcmv1.cpp.

Model expexts input to be processed like cmix -c and some chars swaped. 
By default #TEXTMODE is disabled for this.
If no dictionary is used then #TEXTMODE should be enabled. This is compile time option.

If this model is used then paq8hp should be disabled.
Not tested under Linux so may need adjustments.