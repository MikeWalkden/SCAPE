�$SETUP_MODULE�setup.obj�
3
@(#) FTN95 ver 7.00.0 [MOD ver 7.01] (c) Silverfrost Ltd 1993-2011
SOURCE: C:\workshop\GitHub\SCAPE-MASTER\SCAPE\01 Source code\SCAPE_Inc247\setup.F95
LINE: 40
�
GET_ANGLE_CHANGE K      SETUP_MODULE ��  �       HAVECHANGE    �CHANGE    �YEAR    �
GET_HEIGHT_FACTOR K      SETUP_MODULE ��  �       HAVEFACTOR    �FACTOR    �YEAR    �
GET_PERIOD_FACTOR K      SETUP_MODULE ��  �       HAVEFACTOR    �FACTOR    �YEAR    �
RELATIVE_TO_OD 0      SETUP_MODULE ��  �       ELEMENT    �
SET_OUTPUT_YEAR -      SETUP_MODULE �  �       YEAR    �
DEALLOCATE_PARAMETERS "      SETUP_MODULE �  �        
EQUAL_NOCASE F      SETUP_MODULE ��  �       STRING2  ����  �STRING1  ����  �
UPPERCASE 7      SETUP_MODULE ��  �����       STRING  ����  �
GROYNE_EFFECT_LOOKUP .      SETUP_MODULE �  �       WIDTH    �
OFFSHORE_STR_LOOKUP <      SETUP_MODULE �  �       PERIOD    �HEIGHT    �
CROSSBEACH_STR_LOOKUP <      SETUP_MODULE �  �       PERIOD    �HEIGHT    �
LOG_PARAMETERS "      SETUP_MODULE <  �        
LOG_I_VALUES /      SETUP_MODULE �}  �       VALUES � �
LOG_DP_VALUES /      SETUP_MODULE �{  �       VALUES � �
READ_MODIFY_YEAR_DATA �      SETUP_MODULE .w  �       SUCCESS_CODE    �VALUES  �YEARS  �VALUECOUNT    �KEYWORD  ����  �STEERING_FILE_UNIT    �
INITIALISE_PARAMETERS Q      SETUP_MODULE     �       SUCCESS_CODE    �RADCON    �PATH  ����  �
VELLINGASECTIONSFILE "       SETUP_MODULE �  @     
STARTPROFILEFILE "       SETUP_MODULE    @     
WAVEEXCLUSIONANGLESFILE "       SETUP_MODULE @  @     
GROYNEREMOVALYEARSFILE "       SETUP_MODULE �  @     
GROYNECONSTRUCTIONYEARSFILE "       SETUP_MODULE �  @     
SEAWALLREMOVALYEARSFILE "       SETUP_MODULE    @     
SEAWALLCONSTRUCTIONYEARSFILE "       SETUP_MODULE @  @     
CLIFFDATAFILE "       SETUP_MODULE �  @     
LAYERDATAFILE "       SETUP_MODULE �  @     
BASEGEOLOGYFILE "       SETUP_MODULE    @     
PROFOFFSETSFILE "       SETUP_MODULE @  @     
STARTPROFILESFILE "       SETUP_MODULE �  @     
BEACHCONTENTFILE "       SETUP_MODULE �  @     
ROCKSTRENGTHFILE "       SETUP_MODULE    @     
OFFSHOREBARVOLFILE "       SETUP_MODULE @  @     
TALUSVOLUMEFILE "       SETUP_MODULE �  @     
BEACHVOLUMESFILE "       SETUP_MODULE �  @     
BERMWIDTHSFILE "       SETUP_MODULE    @     
SEALEVELFILE "       SETUP_MODULE @  @     
WAVESANDTIDESFILE "       SETUP_MODULE �  @     
CLIFFELEVATIONSFILE "       SETUP_MODULE �  @     
CERCCONSTANTSFILE "       SETUP_MODULE    @     
OSCANGLESFILE "       SETUP_MODULE @  @     
DZ        SETUP_MODULE `      
DY        SETUP_MODULE h      
DX        SETUP_MODULE p      
OUTPUTCONTOURLEVEL        SETUP_MODULE x      
TALUSSLOPE_RADIANS        SETUP_MODULE �      
BASELINEANGLE_RADIANS        SETUP_MODULE �      
WAVE_EXCLUSION_ANGLE_L_RADIANS        SETUP_MODULE �      
WAVE_EXCLUSION_ANGLE_R_RADIANS        SETUP_MODULE �      
MAXPLATFORMSLOPE_RADIANS        SETUP_MODULE �      
(SETUP_MODULE!MODIFYANGLEAMOUNT:start:1)         �     
(SETUP_MODULE!MODIFYANGLEAMOUNT:size:1)         �     
MODIFYANGLEAMOUNT        SETUP_MODULE     � A
(SETUP_MODULE!MODIFYANGLEYEAR:start:1)         �     
(SETUP_MODULE!MODIFYANGLEYEAR:size:1)         �     
MODIFYANGLEYEAR        SETUP_MODULE    � A
(SETUP_MODULE!MODIFYHEIGHTFACTOR:start:1)         �     
(SETUP_MODULE!MODIFYHEIGHTFACTOR:size:1)         �     
MODIFYHEIGHTFACTOR        SETUP_MODULE    � A
(SETUP_MODULE!MODIFYHEIGHTYEAR:start:1)         �     
(SETUP_MODULE!MODIFYHEIGHTYEAR:size:1)         �     
MODIFYHEIGHTYEAR        SETUP_MODULE    � A
(SETUP_MODULE!MODIFYPERIODFACTOR:start:1)         �     
(SETUP_MODULE!MODIFYPERIODFACTOR:size:1)         �     
MODIFYPERIODFACTOR        SETUP_MODULE    � A
(SETUP_MODULE!MODIFYPERIODYEAR:start:1)         �     
(SETUP_MODULE!MODIFYPERIODYEAR:size:1)         �     
MODIFYPERIODYEAR        SETUP_MODULE    � A
(SETUP_MODULE!GROYNEEFFECTVALUE:start:1)         �     
(SETUP_MODULE!GROYNEEFFECTVALUE:size:1)         �     
GROYNEEFFECTVALUE        SETUP_MODULE    � A
(SETUP_MODULE!GROYNEEFFECTBEACHWIDTH:start:1)         �     
(SETUP_MODULE!GROYNEEFFECTBEACHWIDTH:size:1)         �     
GROYNEEFFECTBEACHWIDTH        SETUP_MODULE     � A
(SETUP_MODULE!OFFSTRVALUE:start:2)         �     
(SETUP_MODULE!OFFSTRVALUE:size:2)         �     
(SETUP_MODULE!OFFSTRVALUE:start:1)         �     
(SETUP_MODULE!OFFSTRVALUE:size:1)         �     
OFFSTRVALUE        SETUP_MODULE (   � B
(SETUP_MODULE!OFFSTRPERIOD:start:2)         �     
(SETUP_MODULE!OFFSTRPERIOD:size:2)         �     
(SETUP_MODULE!OFFSTRPERIOD:start:1)         �     
(SETUP_MODULE!OFFSTRPERIOD:size:1)         �     
OFFSTRPERIOD        SETUP_MODULE 0   � B
(SETUP_MODULE!OFFSTRWAVEHEIGHT:start:1)         �     
(SETUP_MODULE!OFFSTRWAVEHEIGHT:size:1)         �     
OFFSTRWAVEHEIGHT        SETUP_MODULE 8   � A
(SETUP_MODULE!CBSTRVALUE:start:2)         �     
(SETUP_MODULE!CBSTRVALUE:size:2)         �     
(SETUP_MODULE!CBSTRVALUE:start:1)         �     
(SETUP_MODULE!CBSTRVALUE:size:1)         �     
CBSTRVALUE        SETUP_MODULE @   � B
(SETUP_MODULE!CBSTRPERIOD:start:2)         �     
(SETUP_MODULE!CBSTRPERIOD:size:2)         �     
(SETUP_MODULE!CBSTRPERIOD:start:1)               
(SETUP_MODULE!CBSTRPERIOD:size:1)              
CBSTRPERIOD        SETUP_MODULE H   � B
(SETUP_MODULE!CBSTRWAVEHEIGHT:start:1)              
(SETUP_MODULE!CBSTRWAVEHEIGHT:size:1)              
CBSTRWAVEHEIGHT        SETUP_MODULE P   � A
VELLINGAACTIVE        SETUP_MODULE     
DEPTHVALUES        SETUP_MODULE     
PERIODVALUES        SETUP_MODULE     
CELERITYTABLESASTXTFILES        SETUP_MODULE     
INCLUDEGROYNES        SETUP_MODULE      
INCLUDESEAWALL        SETUP_MODULE $    
TIDALRANGEVARIATION        SETUP_MODULE (    
TIDALAMPSCALERIGHT        SETUP_MODULE �     
TIDALAMPSCALELEFT        SETUP_MODULE �     
TIDALAMPSECTION        SETUP_MODULE ,    
MEDIANGRAINSIZE        SETUP_MODULE �     
USEALTSEDIMENTTRANSPORT        SETUP_MODULE 0    
SLUMPPERIOD        SETUP_MODULE 4    
NCLIFFSIMULATIONS        SETUP_MODULE 8    
SWITCHTRANSPORTWHDIFF        SETUP_MODULE <    
WAVEEXCLUSIONANGLER        SETUP_MODULE �     
WAVEEXCLUSIONANGLEL        SETUP_MODULE �     
EXCLUDEWAVES        SETUP_MODULE @    
QPMAXBOUNDARYLEFTOUT        SETUP_MODULE �     
QPMAXBOUNDARYRIGHTOUT        SETUP_MODULE �     
QPMAXBOUNDARYLEFTIN        SETUP_MODULE �     
QPMAXBOUNDARYRIGHTIN        SETUP_MODULE �     
MINBOUNDARYTRANSPORTRATIO        SETUP_MODULE �     
BEACHCRESTLEVEL        SETUP_MODULE �      
BRUUNCONST        SETUP_MODULE      
BEACHHEIGHT        SETUP_MODULE     
MAXBLOCKSIZE        SETUP_MODULE     
TALUSSLOPE        SETUP_MODULE     
TALUSSTRENGTHRATIO        SETUP_MODULE      
BEACHRETURNRATIO        SETUP_MODULE (    
BEACHDISTURBANCERATIO        SETUP_MODULE 0    
MAXBERMSTEP        SETUP_MODULE 8    
MINBERMSTEP        SETUP_MODULE @    
INITIALBERMSTEP        SETUP_MODULE H    
BERMSLOPE        SETUP_MODULE P    
MAXPLATFORMSLOPE        SETUP_MODULE X    
MINHS        SETUP_MODULE `    
MINSTORMWAVEHEIGHT        SETUP_MODULE h    
OUTPUTCONTOURLEVELUSER "       SETUP_MODULE H  
     
DEPTHOSCMSL        SETUP_MODULE p     
MSLOFFSET        SETUP_MODULE x    
MSLM        SETUP_MODULE �    
SEAWALLBASELEVEL        SETUP_MODULE �     
SEAWALLMODE        SETUP_MODULE \    
DEFINESEAWALLSBYYEARANDPOSITION &       SETUP_MODULE `          
DEFINESEAWALLSBYYEAR &       SETUP_MODULE d           
(SETUP_MODULE!DEPTHWPMSL:start:1)         h     
(SETUP_MODULE!DEPTHWPMSL:size:1)         l     
DEPTHWPMSL        SETUP_MODULE X   � A
(SETUP_MODULE!FSWP:start:1)         p     
(SETUP_MODULE!FSWP:size:1)         t     
FSWP        SETUP_MODULE \   � A
NWAVEPOINTS        SETUP_MODULE x    
WAVETIDEFILESTART        SETUP_MODULE |    
LASTACTIVESECTION        SETUP_MODULE �    
FIRSTACTIVESECTION        SETUP_MODULE �    
BASELINEANGLE        SETUP_MODULE �     
NLAYERS        SETUP_MODULE �    
YMAXH        SETUP_MODULE �    
ELEMENTWIDTH        SETUP_MODULE �    
USEVAGRID        SETUP_MODULE �    
ELEMENTHEIGHT        SETUP_MODULE �    
SECTIONHEIGHT        SETUP_MODULE �    
SECTIONWIDTH        SETUP_MODULE �    
NSECTIONS        SETUP_MODULE �    
LASTSECTIONPROFILEOUTPUT        SETUP_MODULE �    
FIRSTSECTIONPROFILEOUTPUT        SETUP_MODULE �    
LASTYEARPROFILEOUTPUT        SETUP_MODULE �    
FIRSTYEARPROFILEOUTPUT        SETUP_MODULE �    
PROFILEOUTPUTTIMESTEP        SETUP_MODULE �    
FIRSTYEARTRANSPORTOUTPUT        SETUP_MODULE �    
HIGHRESOUTSTARTYEAR        SETUP_MODULE �    
HIGHRESOUTTIMESTEP        SETUP_MODULE �    
LOWRESOUTTIMESTEP        SETUP_MODULE �    
BEACHZERODELAY        SETUP_MODULE �    
SPINUPTOYEAR        SETUP_MODULE �    
ENDYEAR        SETUP_MODULE �    
STARTYEAR        SETUP_MODULE �    
RUNMODE        SETUP_MODULE �    
START &       SETUP_MODULE �          
RESTART &       SETUP_MODULE �           
TIDESPERYEAR        SETUP_MODULE �    
DOUBLE         SETUP_MODULE �          