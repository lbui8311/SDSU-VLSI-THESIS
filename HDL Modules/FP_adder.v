`timescale 1ns / 1ps

`define WI_max ((WI1 > WI2) ? WI1 : WI2)
`define WF_max ((WF1 > WF2) ? WF1 : WF2)

module FP_adder #(parameter WI1 = 5, WF1 = 4,        // (5.4)
                  WI2 = 7, WF2 = 3,                  // (7.3) 
                  WIO = 8, WFO = 4)          
                 (input signed [WI1+WF1-1:0] in1,
                  input signed [WI2+WF2-1:0] in2,
                  output signed [WIO+WFO-1:0] out);  // (8.4)
    wire OVF; // Overflow flag
    
    // Temporary placeholders for sign extend and zero pad
    reg signed [(`WI_max+`WF_max):0] in1_temp, in2_temp;
    reg signed [(`WI_max+`WF_max)+1:0] result_temp;
    reg signed [(`WI_max + WIO):0] int_temp;
    reg signed [(`WF_max + WFO):0] frc_temp;

    always @* begin 
        // sign extend WI1 to the length of WI2 & zero pad WF1 to the length of WF2
        in1_temp = { { ((WI2)>(WI1) ? ((WI2)-(WI1)+1) : 1) {in1[(WI1+WF1)-1]}}, in1,
                     { ((WF2)>(WF1) ? ((WF2)-(WF1)) : 0){1'b0}}};
        
        //sign extend WI2 to the length of WI1 & zero pad WF2 to the length of WF1
        in2_temp = { { ((WI1)>(WI2) ? ((WI1)-(WI2)+1) : 1) {in2[(WI2+WF2)-1]}}, in2,
                     { ((WF1)>(WF2) ? ((WF1)-(WF2)) : 0){1'b0}}};   
                     
        // Addition with both inputs sign extended
        result_temp = in1_temp + in2_temp;
        
        // Store WI and WF lengths of result
        int_temp[`WI_max+1:0] = result_temp[(`WI_max+`WF_max)+1-:`WI_max+2];
        frc_temp[`WF_max:0]   = result_temp[`WF_max:0];
        
        // Compare and sign extend resulting WI if less than WIO
        if(WIO > `WI_max+1) int_temp = {{(WIO-`WI_max+1){result_temp[(`WI_max+`WF_max)+1]}}, result_temp[(`WI_max+`WF_max)+1-:`WI_max+2]};
        else                int_temp = {result_temp[(`WI_max+`WF_max)+1], result_temp[(WIO+`WF_max)-1-:WIO]};
        
        // Compare and zero pad resulting WF if less than WFO
        if(WFO > `WF_max) frc_temp = {result_temp[`WF_max:0], {(WFO-`WF_max){1'b0}}};
        else              frc_temp = result_temp[`WF_max-1-:WFO];
        
    end  
       
    // Check to see if dropped bits are the same as new bits
    assign OVF = (result_temp[`WI_max+`WF_max+1:`WI_max+`WF_max] != {(`WI_max-WIO+1){result_temp[`WI_max+`WF_max]}}) ? 1'b0 : 1'b0;
    assign out = {int_temp[WIO+1:0], frc_temp[WFO-1:0]};
endmodule
