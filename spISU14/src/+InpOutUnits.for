************************************************************************
      MODULE InpOutUnits
************************************************************************
*                                                                      *
*     Inputs and Outputs Units                                         *
*                                                                      *
*                                         Angarsk, Russia, 2005/08/22  *
*                                    ITEP, Moscow, Russia, 2011/06/06  *
*                                     JINR, Dubna, Russia, 2012/02/27  *
************************************************************************
*     ---------------------------------------------------------------- *
*     PATHS TO INPUT/OUTPUT FOLDERS                                    *
*     ---------------------------------------------------------------- *
      CHARACTER(*),PARAMETER::
*     GENERAL SETTINGS
     #             Inp  ='../input/',
     #             Out  ='../output/',
*     ATMOSPHERIC NEUTRINO DATA
     #             DatACN=Inp,
     #             OutACN=Out
*     ---------------------------------------------------------------- *
*     UNIT NUMBERS TO INPUT/OUTPUT FILES                               *
*     ---------------------------------------------------------------- *
           INTEGER,PARAMETER::
     #             Nlog  = 10,
     #             Ndata = 11,
     #             Ndat00=201,
     #             Ndat01=Ndat00+01, Ndat02=Ndat00+02,
     #             Ndat03=Ndat00+03, Ndat04=Ndat00+04,
     #             Ndat05=Ndat00+05, Ndat06=Ndat00+06,
     #             Ndat07=Ndat00+07, Ndat08=Ndat00+08,
     #             Ndat09=Ndat00+09, Ndat10=Ndat00+10,
     #             Ndat11=Ndat00+11, Ndat12=Ndat00+12,
     #             Ndat13=Ndat00+13, Ndat14=Ndat00+14,
     #             Ndat15=Ndat00+15, Ndat16=Ndat00+16,
     #             Ndat17=Ndat00+17, Ndat18=Ndat00+18,
     #             Ndat19=Ndat00+19, Ndat20=Ndat00+20,
     #             Ndat21=Ndat00+21, Ndat22=Ndat00+22,
     #             Ndat23=Ndat00+23, Ndat24=Ndat00+24,
     #             Ndat25=Ndat00+25, Ndat26=Ndat00+26,
     #             Ndat27=Ndat00+27, Ndat28=Ndat00+28,
     #             Ndat29=Ndat00+29, Ndat30=Ndat00+30,
     #             Ndat31=Ndat00+31, Ndat32=Ndat00+32,
     #             Ndat33=Ndat00+33, Ndat34=Ndat00+34,
     #             Ndat35=Ndat00+35, Ndat36=Ndat00+36,
     #             Ndat37=Ndat00+37, Ndat38=Ndat00+38,
     #             Ndat39=Ndat00+39, Ndat40=Ndat00+40,
     #             Ndat41=Ndat00+41, Ndat42=Ndat00+42,
     #             Ndat43=Ndat00+43, Ndat44=Ndat00+44,
     #             Ndat45=Ndat00+45, Ndat46=Ndat00+46,
     #             Ndat47=Ndat00+47, Ndat48=Ndat00+48,
     #             Ndat49=Ndat00+49, Ndat50=Ndat00+50,
     #             Ndat51=Ndat00+51, Ndat52=Ndat00+52,
     #             Ndat53=Ndat00+53, Ndat54=Ndat00+54,
     #             Ndat55=Ndat00+55, Ndat56=Ndat00+56,
     #             Ndat57=Ndat00+57, Ndat58=Ndat00+58,
     #             Ndat59=Ndat00+59, Ndat60=Ndat00+60,
     #             Ndat61=Ndat00+61, Ndat62=Ndat00+62,
     #             Ndat63=Ndat00+63, Ndat64=Ndat00+64,
     #             Ndat65=Ndat00+65, Ndat66=Ndat00+66,
     #             Ndat67=Ndat00+67, Ndat68=Ndat00+68,
     #             Ndat69=Ndat00+69, Ndat70=Ndat00+70,
     #             Ndat71=Ndat00+71, Ndat72=Ndat00+72,
     #             Ndat73=Ndat00+73, Ndat74=Ndat00+74,
     #             Ndat75=Ndat00+75, Ndat76=Ndat00+76,
     #             Ndat77=Ndat00+77, Ndat78=Ndat00+78,
     #             Ndat79=Ndat00+79, Ndat80=Ndat00+80,
     #             Ndat81=Ndat00+81, Ndat82=Ndat00+82,
     #             Ndat83=Ndat00+83, Ndat84=Ndat00+84,
     #             Ndat85=Ndat00+85, Ndat86=Ndat00+86,
     #             Ndat87=Ndat00+87, Ndat88=Ndat00+88,
     #             Ndat89=Ndat00+89, Ndat90=Ndat00+90,
     #             Ndat91=Ndat00+91, Ndat92=Ndat00+92,
     #             Ndat93=Ndat00+93, Ndat94=Ndat00+94,
     #             Ndat95=Ndat00+95, Ndat96=Ndat00+96,
     #             Ndat97=Ndat00+97, Ndat98=Ndat00+98,
     #             Ndat99=Ndat00+99
*     ---------------------------------------------------------------- *

      END MODULE InpOutUnits
