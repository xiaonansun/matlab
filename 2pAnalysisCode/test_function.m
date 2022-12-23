function test_function
disp(mfilename)

st = dbstack;
namestr = st.name;

disp(['DBstack name is: ' namestr])
