�$SLR_MODULE�SeaLevelRise.obj�
3
@(#) FTN95 ver 7.00.0 [MOD ver 7.01] (c) Silverfrost Ltd 1993-2011
SOURCE: C:\workshop\GitHub\SCAPE-MASTER\SCAPE\01 Source code\SCAPE_Inc247\SeaLevelRise.F95
LINE: 39
�
SL_FOR_YEAR G      SLR_MODULE     �       SUCCESS_CODE    �SL    �YEAR    �
SLYEAR        SLR_MODULE        
SLFILEUNIT        SLR_MODULE      
SEALEVELFILE "       SETUP_MODULE @  @      
SUCCESS_CODE_FROM_IOSTAT J      EXCEPTIONS     �       MESSAGE  ����  �IOSTAT    �UNIT    �
SUCCESS_CODE_FROM_ALLOCATION D      EXCEPTIONS     �       MESSAGE  ����  �ALLOC_ERROR    �
RAISE_EXCEPTION 2      EXCEPTIONS     �       MESSAGE  ����  �
FELOGSTRING         LOCAL_DATA 4        
FILE_PATH         LOCAL_DATA 4         
BINARYOPENMODE "       UTILITIES_MODULE       
ACCESSMODE &       UTILITIES_MODULE          
FORMMODE &       UTILITIES_MODULE           
DO_OPEN_TEXT_OUTPUT_FILE y      UTILITIES_MODULE     �       UNITS  ����H  �VALUESPERYEAR  H  �FILENAME  ����  �PATH  ����  �UNITNO    �
DO_OPEN_TEXT_INPUT_FILE U      UTILITIES_MODULE     �       FILENAME  ����  �PATH  ����  �UNITNO    �
DO_OPEN_BINARY_INPUT_FILE U      UTILITIES_MODULE     �       FILENAME  ����  �PATH  ����  �UNITNO    �
OPEN_BINARY_FILE P      UTILITIES_MODULE     �       IOS    �FILENAME  ����  �UNITNO    �
ELAPSED_TIME x      UTILITIES_MODULE     �       SECONDS    �MINUTES    �HOURS    �DAYS    �STOP_DAY    �START_DAY    �
JULIAN_DAY l      UTILITIES_MODULE     �       SECOND    �MINUTE    �HOUR    �YEAR    �DAY    �MONTH    �
GREGORIAN_DATE X      UTILITIES_MODULE     �       YEAR    �DAY    �MONTH    �JULIAN_DAY    �
TO_MJD /      UTILITIES_MODULE     �       JD    �
TO_JD 0      UTILITIES_MODULE     �       MJD    �
GET_FILE_UNIT &      UTILITIES_MODULE     �        
SCAPE_LOG 5      UTILITIES_MODULE     �       TOSCREEN  @  �
FELOG &      UTILITIES_MODULE     �        
INTERPOLATE_VALUES K      UTILITIES_MODULE     �       VALUE2     �VALUE1     �DEST � �
AVERAGE_VALUE 3      UTILITIES_MODULE     �       VALUES � �
MASKED_AVERAGE_VALUE        UTILITIES_MODULE     �       MASK � �VALUES � �