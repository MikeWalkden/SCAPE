�$VERTICAL_GRID�vertical_grid.obj�
3
@(#) FTN95 ver 7.00.0 [MOD ver 7.01] (c) Silverfrost Ltd 1993-2011
SOURCE: C:\workshop\GitHub\SCAPE-MASTER\SCAPE\01 Source code\SCAPE_Inc247\vertical_grid.F95
LINE: 51
�
FINISH_VA_GRID #      VERTICAL_GRID /@  �        
SAVE_PROFILE_POINTS 6      VERTICAL_GRID �=  �       SUCCESS_CODE    �
ADD_PROFILE_POINT F      VERTICAL_GRID 59  �       SUCCESS_CODE    �Z    �Y    �
INIT_PROFILE_POINT_STORAGE 6      VERTICAL_GRID 58  �       SUCCESS_CODE    �
LAYER_PROFILE_VALUE J      VERTICAL_GRID �6  �       SECTION    �LAYER    �YINDEX    �
BASE_PROFILE_VALUE >      VERTICAL_GRID �5  �       SECTION    �YINDEX    �
GET_CELL_DOMINANT_LAYER K      VERTICAL_GRID �0  �       SECTION    �YINDEX    �ZVALUE    �
GET_CELL_VALUES i      VERTICAL_GRID �/  �       SEDIMENT    �STRENGTH    �SECTION    �YINDEX    �ZVALUE    �
ADD_MATERIAL_TO_LAYER ~      VERTICAL_GRID T,  �       SUCCESS_CODE    �YPOSITION    �THICKNESS    �VOLUME    �SECTION    �NLAYER    �
ADD_NEW_LAYER �      VERTICAL_GRID �*  �       SUCCESS_CODE    �YPOSITION    �THICKNESS    �VOLUME    �SEDIMENT    �STRENGTH    �SECTION    �
CUT_LAYER_PROFILES K      VERTICAL_GRID �&  �       CLIFFY  �NCE    �NYSECTIONS    �
GET_Y_VALUES_FROM_PROFILES K      VERTICAL_GRID �!  �       CLIFFY  �NCE    �NYSECTIONS    �
GET_BASE_PROFILE 6      VERTICAL_GRID �  �       SUCCESS_CODE    �
INTERPOLATE_Z_VALUE O      VERTICAL_GRID Z  �       Z2    �Y2    �Z1    �Y1    �Y    �
INTERPOLATE_Y_VALUE O      VERTICAL_GRID   �       Z2    �Y2    �Z1    �Y1    �Z    �
INTERPOLATE_PROFILE_VALUE a      VERTICAL_GRID �  �       Z2    �Y2    �Z1    �Y1    �INDEX    �PROFILE � �
INTERPOLATE_PROFILE j      VERTICAL_GRID   �       SUCCESS_CODE    �PROFILE � �PROFILE_END    �PROFILE_START    �
READ_PROFILE �      VERTICAL_GRID   �       SUCCESS_CODE    �BASE_GEOLOGY    �PROFILE_END    �PROFILE_START2    �PROFILE_START    �FILENAME  ����  �UNITNO    �
GET_LAYERS 6      VERTICAL_GRID     �       SUCCESS_CODE    �
NEXT_INDEX        VERTICAL_GRID       
INDEX_LIMIT        VERTICAL_GRID        
ACTIVE_TEMP_POINTS        VERTICAL_GRID $      
(VERTICAL_GRID!TEMP_Z_2:start:1)         (      
(VERTICAL_GRID!TEMP_Z_2:size:1)         ,      
TEMP_Z_2         VERTICAL_GRID     � A
(VERTICAL_GRID!TEMP_Y_2:start:1)         0      
(VERTICAL_GRID!TEMP_Y_2:size:1)         4      
TEMP_Y_2         VERTICAL_GRID    � A
(VERTICAL_GRID!TEMP_Z_1:start:1)         8      
(VERTICAL_GRID!TEMP_Z_1:size:1)         <      
TEMP_Z_1         VERTICAL_GRID    � A
(VERTICAL_GRID!TEMP_Y_1:start:1)         @      
(VERTICAL_GRID!TEMP_Y_1:size:1)         D      
TEMP_Y_1         VERTICAL_GRID    � A
NHCELLS        VERTICAL_GRID H      
DOUBLE &       SETUP_MODULE �          
NLAYERS        SETUP_MODULE �     
DY        SETUP_MODULE h      
DZ         SETUP_MODULE `      