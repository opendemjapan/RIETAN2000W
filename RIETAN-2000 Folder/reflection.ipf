Macro reflection()
	if (strsearch(CsrWave(A),"yphase_1",0) == 0)
		Print "hkl  =",h_1(pcsr(A)),k_1(pcsr(A)),l_1(pcsr(A))," d  =",d_1(pcsr(A))," 2-theta  =",xphase_1(pcsr(A))
	endif
	if (strsearch(CsrWave(A),"yphase_2",0) == 0)
		Print "hkl  =",h_2(pcsr(A)),k_2(pcsr(A)),l_2(pcsr(A))," d  =",d_2(pcsr(A))," 2-theta  =",xphase_2(pcsr(A))
	endif
	if (strsearch(CsrWave(A),"yphase_3",0) == 0)
		Print "hkl  =",h_3(pcsr(A)),k_3(pcsr(A)),l_3(pcsr(A))," d  =",d_3(pcsr(A))," 2-theta  =",xphase_3(pcsr(A))
	endif
	if (strsearch(CsrWave(A),"yphase_4",0) == 0)
		Print "hkl  =",h_4(pcsr(A)),k_4(pcsr(A)),l_4(pcsr(A))," d  =",d_4(pcsr(A))," 2-theta  =",xphase_4(pcsr(A))
	endif
	if (strsearch(CsrWave(A),"yphase_5",0) == 0)
		Print "hkl  =",h_5(pcsr(A)),k_5(pcsr(A)),l_5(pcsr(A))," d  =",d_5(pcsr(A))," 2-theta  =",xphase_5(pcsr(A))
	endif
	if (strsearch(CsrWave(A),"yphase_6",0) == 0)
		Print "hkl  =",h_6(pcsr(A)),k_6(pcsr(A)),l_6(pcsr(A))," d  =",d_6(pcsr(A))," 2-theta  =",xphase_6(pcsr(A))
	endif
	if (strsearch(CsrWave(A),"yphase_7",0) == 0)
		Print "hkl  =",h_7(pcsr(A)),k_7(pcsr(A)),l_7(pcsr(A))," d  =",d_7(pcsr(A))," 2-theta  =",xphase_7(pcsr(A))
	endif
	if (strsearch(CsrWave(A),"yphase_8",0) == 0)
		Print "hkl  =",h_8(pcsr(A)),k_8(pcsr(A)),l_8(pcsr(A))," d  =",d_8(pcsr(A))," 2-theta  =",xphase_8(pcsr(A))
	endif
end
