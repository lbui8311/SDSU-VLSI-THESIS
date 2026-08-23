`timescale 1ns / 1ps


module FP_Inv_Mitchell
#(
    parameter WI = 5, WF = 11,
    parameter WIO = 2, WFO = 14
)(
    input RST,
    input [(WI+WF)-1:0] x,                // Q(5.11)
    output reg signed [(WIO+WFO)-1:0] y       // Q(2.14)
); 

    localparam ONE = 1 << WF;
    
    integer e;
    
    reg [WF-1:0] F;
    reg [2* (WI+WF)-1:0] M, temp;
    
    always @* begin
        e = x >>> WF;
        F = x[WF-1:0];
        M = ONE | M;
        
        if (e >= (WI+WF))       temp = {2*(WI+WF){1'b1}};
        else if (e < -(WI+WF))  temp = {2*(WI+WF){1'b0}};
        else if (e >= 0)        temp = M << e;
        else                    temp = M >> (-e);
        
        y = temp[(WI+WF)-1:0];
    end
    
endmodule
