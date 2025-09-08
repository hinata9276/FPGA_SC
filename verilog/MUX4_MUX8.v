module mux4 (
    input [3:0] data,  // 4-bit input data (D0, D1, D2, D3)
    input [1:0] select, // 2-bit select input (S1, S0)
    output reg out     // Output of the multiplexer
);

    always @(*) begin
        case (select)
            2'b00: out = data[0]; // Select D0
            2'b01: out = data[1]; // Select D1
            2'b10: out = data[2]; // Select D2
            2'b11: out = data[3]; // Select D3
            default: out = 1'b0;  // Default value (optional)
        endcase
    end

endmodule

module mux8 (
    input [7:0] data,    // 8-bit input data (D0 to D7)
    input [2:0] select,  // 3-bit select input (S2, S1, S0)
    output out           // Output of the multiplexer
);
    wire mux4_out0;  // Output from the first MUX4
    wire mux4_out1;  // Output from the second MUX4

    // Instantiate two MUX4 modules
    mux4 mux4_0 (
        .data(data[3:0]),         // Lower 4 inputs (D0 to D3)
        .select(select[1:0]),     // Least significant select bits
        .out(mux4_out0)
    );

    mux4 mux4_1 (
        .data(data[7:4]),         // Upper 4 inputs (D4 to D7)
        .select(select[1:0]),     // Least significant select bits
        .out(mux4_out1)
    );

    MUXF7 MUXF7_inst (
        .O(out),    // Output of MUX to general routing
        .I0(mux4_out0),  // Input (tie to LUT6 O6 pin)
        .I1(mux4_out1),  // Input (tie to LUT6 O6 pin)
        .S(select[2])     // Input select to MUX
    );

endmodule