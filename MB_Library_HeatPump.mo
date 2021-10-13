package MB_Library_HeatPump

  package Heat_Pump
    model ExpanValve
      Modelica.Blocks.Interfaces.RealInput Pc
        "Pressure in the condensor"
        annotation (Placement(transformation(extent={{-120,-30},{-100,-10}})));
      Modelica.Blocks.Interfaces.RealInput Pe
        "Pressure in the evaporator"
        annotation (Placement(transformation(extent={{-120,-70},{-100,-50}})));
      Modelica.Blocks.Interfaces.RealOutput MassFlow_out
        annotation (Placement(transformation(extent={{100,50},{120,70}})));
      Modelica.Blocks.Interfaces.RealOutput H_out
        "Specific enthalpy output of the expansion valve"
        annotation (Placement(transformation(extent={{100,10},{120,30}})));
      Modelica.Blocks.Interfaces.RealInput H_in
        "Specific enthalpy input of the expansion valve"
        annotation (Placement(transformation(extent={{-120,50},{-100,70}})));
      Modelica.Blocks.Interfaces.RealInput A_v "Opening area" annotation (
          Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=90,
            origin={0,-110})));
      Modelica.Blocks.Interfaces.RealInput D_in "Density of the saturated liquid"
        annotation (Placement(transformation(extent={{-120,10},{-100,30}})));

    Real Cv;

    equation

    MassFlow_out = A_v * Cv * sqrt(D_in*(Pc - Pe));
    Cv = 3 + 5e-6 * (Pc - Pe)                           "Orifice coefficient";
    H_in = H_out                                        "Adiabatic";

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              lineThickness=1), Bitmap(
              extent={{-92,-90},{94,94}},
              imageSource="PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iaXNvLTg4NTktMSI/Pg0KPCEtLSBHZW5lcmF0b3I6IEFkb2JlIElsbHVzdHJhdG9yIDE5LjAuMCwgU1ZHIEV4cG9ydCBQbHVnLUluIC4gU1ZHIFZlcnNpb246IDYuMDAgQnVpbGQgMCkgIC0tPg0KPHN2ZyB2ZXJzaW9uPSIxLjEiIGlkPSJDYXBhXzEiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgeG1sbnM6eGxpbms9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGxpbmsiIHg9IjBweCIgeT0iMHB4Ig0KCSB2aWV3Qm94PSIwIDAgNDgwIDQ4MCIgc3R5bGU9ImVuYWJsZS1iYWNrZ3JvdW5kOm5ldyAwIDAgNDgwIDQ4MDsiIHhtbDpzcGFjZT0icHJlc2VydmUiPg0KPGc+DQoJPGc+DQoJCTxwYXRoIGQ9Ik00MTYsMjY3Ljk5OHYtMTZoLTMydi0yNGgtNDh2NDBoLTI4LjE2OGMtNC4zOTEtNS4zODYtOS40NDItMTAuMi0xNS4wMzItMTQuMzI4Yy0xLjYtMS4yLTMuMi0yLjMyLTQuOC0zLjM1MnYtMzguMzJoMjQNCgkJCXYtNDhoLTQwdi0xNmgyNHYtMzJoNjRjNC40MTgsMCw4LTMuNTgyLDgtOHYtMTZjMC0xMy4yNTUtMTAuNzQ1LTI0LTI0LTI0SDEzNmMtMTMuMjU1LDAtMjQsMTAuNzQ1LTI0LDI0djE2YzAsNC40MTgsMy41ODIsOCw4LDgNCgkJCWg2NHYzMmgyNHYxNmgtNDB2NDhoMjR2MzguMzJjLTEuNiwxLjAzMi0zLjIsMi4xMzYtNC43NDQsMy4zMmMtNS42MTQsNC4xMzItMTAuNjgzLDguOTU3LTE1LjA4OCwxNC4zNkgxNDR2LTQwSDk2djI0SDY0djE2SDB2MTEyDQoJCQloNjR2MTZoMzJ2MTZoNDh2LTMyaDI4LjE4NGMzMC44NDcsMzcuNDU0LDg2LjIxNiw0Mi44MDksMTIzLjY3LDExLjk2MmM0LjM2My0zLjU5Myw4LjM2OS03LjU5OSwxMS45NjItMTEuOTYySDMzNnYzMmg0OHYtMTZoMzINCgkJCXYtMTZoNjR2LTExMkg0MTZ6IE02NCwzNjMuOTk4SDE2di04MGg0OFYzNjMuOTk4eiBNOTYsMzc5Ljk5OEg4MHYtMTEyaDE2VjM3OS45OTh6IE0xMjgsMzk1Ljk5OGgtMTZ2LTE1MmgxNlYzOTUuOTk4eg0KCQkJIE0xMjgsOTkuOTk4di04YzAtNC40MTgsMy41ODItOCw4LThoMjA4YzQuNDE4LDAsOCwzLjU4Miw4LDh2OEgxMjh6IE0yMDAsMTMxLjk5OHYtMTZoODB2MTZIMjAweiBNMjU2LDE0Ny45OTh2MTZoLTMydi0xNkgyNTZ6DQoJCQkgTTE4NCwxOTUuOTk4di0xNmgxMTJ2MTZIMTg0eiBNMzM2LDM2My45OThoLTMyLjA4Yy0yLjUxOCwwLTQuODg5LDEuMTg2LTYuNCwzLjJjLTIzLjc3LDMxLjc2Ny02OC43OTMsMzguMjUtMTAwLjU2LDE0LjQ4DQoJCQljLTUuNDkzLTQuMTExLTEwLjM3LTguOTg3LTE0LjQ4LTE0LjQ4Yy0xLjUxMS0yLjAxNC0zLjg4Mi0zLjItNi40LTMuMkgxNDR2LTgwaDMyLjA4YzIuNTI4LTAuMDEzLDQuOTAxLTEuMjIsNi40LTMuMjU2DQoJCQljNS44MTYtNy43MiwxMy4xMzctMTQuMTgxLDIxLjUyLTE4Ljk5MmMyLjQ4My0xLjQzMyw0LjAwOS00LjA4NSw0LTYuOTUydi00Mi44aDY0djQyLjhjLTAuMDA4LDIuODY3LDEuNTE3LDUuNTE5LDQsNi45NTINCgkJCWMyLjUsMS40NDYsNC45MDQsMy4wNDksNy4yLDQuOGM1LjQzLDMuOTk4LDEwLjIyNiw4Ljc5NCwxNC4yMjQsMTQuMjI0YzEuNDk5LDIuMDM2LDMuODcyLDMuMjQzLDYuNCwzLjI1NkgzMzZWMzYzLjk5OHoNCgkJCSBNMzY4LDM5NS45OThoLTE2di0xNTJoMTZWMzk1Ljk5OHogTTQwMCwzNzkuOTk4aC0xNnYtMTEyaDE2VjM3OS45OTh6IE00NjQsMzYzLjk5OGgtNDh2LTgwaDQ4VjM2My45OTh6Ii8+DQoJPC9nPg0KPC9nPg0KPGc+DQo8L2c+DQo8Zz4NCjwvZz4NCjxnPg0KPC9nPg0KPGc+DQo8L2c+DQo8Zz4NCjwvZz4NCjxnPg0KPC9nPg0KPGc+DQo8L2c+DQo8Zz4NCjwvZz4NCjxnPg0KPC9nPg0KPGc+DQo8L2c+DQo8Zz4NCjwvZz4NCjxnPg0KPC9nPg0KPGc+DQo8L2c+DQo8Zz4NCjwvZz4NCjxnPg0KPC9nPg0KPC9zdmc+DQo=",
              fileName="//nas.ads.mwn.de/ge43kok/Desktop/SVGs/valve.svg")}),
                                                                     Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end ExpanValve;

    model Condensor_R134a

      Modelica.Blocks.Interfaces.RealInput H_refri_in        "Inlet refrigerant specific enthalpy" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-110,0})));
      Modelica.Blocks.Interfaces.RealInput MassFlow_refri_in       "Inlet refrigerant massflow"   annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-110,40})));
      Modelica.Blocks.Interfaces.RealInput MassFlow_refri_out      "Outlet refrigerant massflow"   annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-110,80})));
      Modelica.Blocks.Interfaces.RealOutput H_refri_out      "Outlet refrigerant specific enthalpy" annotation (Placement(
            transformation(extent={{100,50},{120,70}}), iconTransformation(extent=
               {{100,50},{120,70}})));
      Modelica.Blocks.Interfaces.RealInput T_water_in        "Inlet condensing water temperature" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-110,-40})));
      Modelica.Blocks.Interfaces.RealOutput Pc "Refrigerant pressure in the condensor"                 annotation (Placement(transformation(extent={{100,-30},{120,-10}}),
            iconTransformation(extent={{100,-30},{120,-10}})));
      Modelica.Blocks.Interfaces.RealInput MassFlow_water_in "Inlet condensing water massflow" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-110,-80})));
      Modelica.Blocks.Interfaces.RealOutput Power_heat       "Heat power" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={40,110}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={40,110})));
      Modelica.Blocks.Interfaces.RealOutput T_water_out      "Outlet condensing water temperauture" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-40,110})));
      Modelica.Blocks.Interfaces.RealOutput D_refri_out      "Outlet refrigerant density" annotation (Placement(transformation(extent={{100,10},{120,30}}),
            iconTransformation(extent={{100,10},{120,30}})));

    import si = Modelica.SIunits;

    /************************* Define all parameters ***************************************/
            parameter Real n = 1                                       "Number of tubes"                                        annotation (Dialog(tab="Physical parameters"));
            parameter si.Length Ro = 12.5e-3                           "Outer radius of the tube"                               annotation (Dialog(tab="Physical parameters"));
            parameter si.Length Ri = 11.5e-3                           "Inner radius of the tube"                               annotation (Dialog(tab="Physical parameters"));
            parameter si.Length Lt = 15                                "Total length of each tube"                              annotation (Dialog(tab="Physical parameters"));
            parameter si.CoefficientOfHeatTransfer Alpa_i_vapor = 400  "Inner heat transfer coefficient of the vapor zone"      annotation (Dialog(tab="Physical parameters"));
            parameter si.CoefficientOfHeatTransfer Alpa_i_satu = 2000  "Inner heat transfer coefficient of the two phase zone"  annotation (Dialog(tab="Physical parameters"));
            parameter si.CoefficientOfHeatTransfer Alpa_o = 2000       "Outer heat transfer coefficient of tubes for all zones" annotation (Dialog(tab="Physical parameters"));
            parameter si.SpecificHeatCapacity Cp_wall = 385            "Specific heat of pipe wall material"                    annotation (Dialog(tab="Physical parameters"));
            parameter si.Density D_wall = 8960                         "Density of pipe wall material"                          annotation (Dialog(tab="Physical parameters"));
            parameter si.SpecificHeatCapacity  Cp_water = 4185         "Specific heat capacity of water"                        annotation (Dialog(tab="Physical parameters"));
            parameter si.Density D_water = 993                         "Density of water"                                       annotation (Dialog(tab="Physical parameters"));
            parameter si.Volume Volume_water = 0.37                    "The volume for condensing water in the condensor"       annotation (Dialog(tab="Physical parameters"));
            parameter si.Pressure Pc_initial = 6e5                     "Initial refrigerant pressure in the condensor" annotation (Dialog(tab="Initial conditions"));
            parameter si.SpecificEnthalpy H_refri_out_initial = 320e3  "Initial output refrigerant specific enthalpy of the condensor"   annotation (Dialog(tab="Initial conditions"));
            parameter si.Temperature T_wall_1_initial = 293.15         "For condensor, initial wall temperature of the vapor zone"       annotation (Dialog(tab="Initial conditions"));
            parameter si.Temperature T_wall_2_initial = 293.15         "For condensor, initial wall temperature of the two-phase zone"   annotation (Dialog(tab="Initial conditions"));

     /************************* Define all variables ***************************************/
            Real Gama                         "Average void fraction";
            Real a;

            si.Length L1                      "For condensor, length of the vapor zone";
            si.Temperature T_wall_1           "For condensor, wall temperature of the vapor zone";
            si.Temperature T_wall_2           "For condensor, wall temperature of the two-phase zone";
            si.Temperature T_wall_boundary    "Temperature at the boundary between two zones";
            si.MassFlowRate MassFlow_boundary "Massflow at the boundary between two zones";

            Modelica.Media.R134a.R134a_ph.ThermodynamicState Vapor_state  "Thermal dynamic state variable for refrigerant vapor";
            Modelica.Media.R134a.R134a_ph.SaturationProperties Satu_state "Thermal dynamic state variable for two-phase refrigerant mixture";

            si.Area A_cs_o "Tube outer crossection area (pipe thickness counted)";
            si.Area A_cs_i "Tube inner crossection area (hydraulic)";
            si.Area A_sf_i "Tube inner surface area (contact area between refrigerant and pipe wall)";
            si.Area A_sf_o "Tube outer surface area (contact area between air and pipe wall)";

            si.Temperature T_refri_satu             "Refrigerant temperature in two-phase zone";
            si.SpecificEnthalpy H_refri_satu_liquid "Refrigerant specific enthalpy for saturated liquid state";
            si.SpecificEnthalpy H_refri_satu_gas    "Refrigerant specific enthalpy for saturated gas state";
            si.Density D_refri_satu_liquid          "Refrigerant density for saturated liquid state";
            si.Density D_refri_satu_gas             "Refrigerant density for saturated gas state";

            si.SpecificEnthalpy H_refri_vapor_avg   "Refrigerant specific enthalpy for vapor zone (averaged)";
            si.Density D_refri_vapor_avg            "Refrigerant density for vapor zone (averaged)";
            si.Temperature T_refri_vapor_avg        "Refrigerant temperature for vapor zone (averaged)";

            Real dDdH_refri_at_Pc_vapor                   "Refrigerant thermal dynamics property, calculated from data";
            Real dDdP_refri_at_H_vapor                    "Refrigerant thermal dynamics property, calculated from data";
            Real dDdP_refri_satu_liquid                   "Refrigerant thermal dynamics property, calculated from data";
            Real dDdP_refri_satu_gas                      "Refrigerant thermal dynamics property, calculated from data";
            Real dHdP_refri_satu_liquid                   "Refrigerant thermal dynamics property, calculated from data";
            Real dHdP_refri_satu_gas                      "Refrigerant thermal dynamics property, calculated from data";

    initial equation

            Pc = Pc_initial;
            H_refri_out = H_refri_out_initial;
            T_wall_1 = T_wall_1_initial;
            T_wall_2 = T_wall_2_initial;
            T_water_out = T_water_in;

    equation
            A_cs_o = Ro^2 * Modelica.Constants.pi "Outer cross-section area of condensor tubes";
            A_cs_i = Ri^2 * Modelica.Constants.pi "Inner cross-section area of condensor tubes";
            A_sf_i = n * Ri * 2 * Modelica.Constants.pi * Lt "Inner heat transfer effective area ";
            A_sf_o = n * Ro * 2 * Modelica.Constants.pi * Lt "Outer heat transfer effective area ";

            a = (D_refri_satu_gas/D_refri_satu_liquid)^(2/3);
            Gama = (( 1 - a * (1 - log(a)))   / (1-a)^2);

    /************************* Fluid properties calculation for the two-phase zone ***************/
            T_refri_satu = Modelica.Media.R134a.R134a_ph.saturationTemperature(Pc);
            Satu_state.psat = Pc;
            Satu_state.Tsat = T_refri_satu;
            dDdP_refri_satu_liquid  = Modelica.Media.R134a.R134a_ph.dBubbleDensity_dPressure(Satu_state);
            dDdP_refri_satu_gas = Modelica.Media.R134a.R134a_ph.dDewDensity_dPressure(Satu_state);
            dHdP_refri_satu_liquid = Modelica.Media.R134a.R134a_ph.dBubbleEnthalpy_dPressure(Satu_state);
            dHdP_refri_satu_gas = Modelica.Media.R134a.R134a_ph.dDewEnthalpy_dPressure(Satu_state);
            H_refri_satu_liquid = Modelica.Media.R134a.R134a_ph.bubbleEnthalpy(Satu_state);
            H_refri_satu_gas = Modelica.Media.R134a.R134a_ph.dewEnthalpy(Satu_state);
            D_refri_satu_liquid = Modelica.Media.R134a.R134a_ph.bubbleDensity(Satu_state);
            D_refri_satu_gas  = Modelica.Media.R134a.R134a_ph.dewDensity(Satu_state);

    /************************* Fluid properties calculation for the vapor zone ***************/
            Vapor_state.phase = 1;
            Vapor_state.h = H_refri_vapor_avg;
            Vapor_state.d = D_refri_vapor_avg;
            Vapor_state.T = T_refri_vapor_avg;
            Vapor_state.p = Pc;
            D_refri_vapor_avg = Modelica.Media.R134a.R134a_ph.density_ph(Pc, H_refri_vapor_avg, 1);
            T_refri_vapor_avg = Modelica.Media.R134a.R134a_ph.temperature_ph(Pc, H_refri_vapor_avg, 1);
            dDdH_refri_at_Pc_vapor  = Modelica.Media.R134a.R134a_ph.density_derh_p(Vapor_state);
            dDdP_refri_at_H_vapor =  Modelica.Media.R134a.R134a_ph.density_derp_h(Vapor_state);

            T_wall_boundary = T_wall_1                                                                                               "The boundary temperature between two-phase and vapor zone is set to be the two-phase zone temperature";
            H_refri_vapor_avg = 0.5 * (H_refri_satu_gas + H_refri_in)                                                                "The refrigerant specific enthalpy in vapor zone is an averaged valve between input and output to the zone";
            D_refri_out = Gama  * D_refri_satu_gas  + (1 - Gama)  * D_refri_satu_liquid                                              "The output refrigerant density is the averaged density in the two-phase zone";
            Power_heat = Alpa_o*A_sf_o*(L1/Lt)*(T_wall_1 - T_water_out) + Alpa_o*A_sf_o*((Lt - L1)/Lt)*(T_wall_2 - T_water_out)      "Heat power is the air-to-wall heat exchange";

    /************************* "Conservation of mass in the vapor region" *****************************************/
           (dDdP_refri_at_H_vapor + 0.5 * dDdH_refri_at_Pc_vapor * dHdP_refri_satu_gas) * A_cs_i * L1 * der(Pc)
           + 0.5* dDdH_refri_at_Pc_vapor * A_cs_i * L1 * der(H_refri_in)
           + (D_refri_vapor_avg - D_refri_satu_gas) * A_cs_i * der(L1)
           =
           MassFlow_refri_in
           - MassFlow_boundary;

    /************************* "Conservation of mass in two-phase region" *****************************************/
           (dDdP_refri_satu_liquid  * (1 - Gama)  + dDdP_refri_satu_gas * Gama)  * A_cs_i * (Lt - L1) * der(Pc)
           + (D_refri_satu_gas  - D_refri_satu_liquid) * (1 - Gama)  * A_cs_i * der(L1)
           =
           MassFlow_boundary
           - MassFlow_refri_out;
           //+ (D_refri_satu_liquid - D_refri_satu_gas) * A_cs_i * (Lt - L1) * der(Gama)

    /************************* "Conservation of energy in vapor region" *******************************************/
           ((dDdP_refri_at_H_vapor + 0.5*dHdP_refri_satu_gas*dDdH_refri_at_Pc_vapor) * H_refri_vapor_avg + 0.5* dHdP_refri_satu_gas * D_refri_vapor_avg - 1) * A_cs_i * L1 * der(Pc)
           + 0.5 * (dDdH_refri_at_Pc_vapor * H_refri_vapor_avg + D_refri_vapor_avg) * A_cs_i * L1 * der(H_refri_in)
           + (D_refri_vapor_avg * H_refri_vapor_avg - D_refri_satu_gas * H_refri_satu_gas) * A_cs_i * der(L1)
           =
           MassFlow_refri_in * H_refri_in
           - MassFlow_boundary * H_refri_satu_gas
           + Alpa_i_vapor * A_sf_i * (L1/Lt) * (T_wall_1 - T_refri_vapor_avg);

    /************************* "Conservation of energy in two-phase region"****************************************/
          ((dDdP_refri_satu_liquid  * H_refri_satu_liquid + dHdP_refri_satu_liquid * D_refri_satu_liquid) * (1 - Gama)  + (dDdP_refri_satu_gas * H_refri_satu_gas + dHdP_refri_satu_gas * D_refri_satu_gas) * Gama - 1) * A_cs_i * (Lt - L1) * der(Pc)
          + (1 - Gama)  * (D_refri_satu_gas  * H_refri_satu_gas - D_refri_satu_liquid * H_refri_satu_liquid) * A_cs_i * der(L1)
          =
          MassFlow_boundary * H_refri_satu_gas
          - MassFlow_refri_out * H_refri_out
          + Alpa_i_satu * A_sf_i *((Lt - L1)/Lt) * (T_wall_2 - T_refri_satu);
          //+ (D_refri_satu_liquid * H_refri_satu_liquid - D_refri_satu_gas * H_refri_satu_gas) * A_cs_i * (Lt - L1) * der(Gama)

    /************************* "Conservation of energy of pipe wall (vapor region)"********************************/
          Cp_wall * D_wall * (A_cs_o - A_cs_i) * L1 * der(T_wall_1)
          =
          Alpa_i_vapor * A_sf_i * (L1/Lt) * (T_refri_vapor_avg - T_wall_1)
          + Alpa_o * A_sf_o * (L1/Lt) * (T_water_out - T_wall_1);

    /************************* "Conservation of energy of pipe wall (two-phase region)"****************************/
          Cp_wall * D_wall * (A_cs_o - A_cs_i) * (Lt - L1) * (der(T_wall_2) - ((T_wall_2 - T_wall_boundary)/(Lt - L1)) * der(L1))
          =
          Alpa_i_satu * A_sf_i * ((Lt - L1)/Lt) * (T_refri_satu - T_wall_2)
          + Alpa_o * A_sf_o * ((Lt - L1)/Lt) * (T_water_out - T_wall_2);

    /************************* "Conservation of energy of condensing water volume "********************************/
          Cp_water * Volume_water * D_water * der(T_water_out)
          =
          Alpa_o * A_sf_o * (L1/Lt) * (T_wall_1 - T_water_out)
          + Alpa_o * A_sf_o * ((Lt - L1)/Lt) * (T_wall_2 - T_water_out)
          + MassFlow_water_in * Cp_water * (T_water_in - T_water_out);

        annotation (Placement(transformation(extent={{100,-70},{120,-50}})),
                  Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Bitmap(
              extent={{-100,-100},{98,98}},
              imageSource="PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiIHN0YW5kYWxvbmU9Im5vIj8+CjwhLS0gQ3JlYXRlZCB3aXRoIElua3NjYXBlIChodHRwOi8vd3d3Lmlua3NjYXBlLm9yZy8pIC0tPgoKPHN2ZwogICB3aWR0aD0iNTIuNDY2ODYybW0iCiAgIGhlaWdodD0iNjEuNjM2MTkybW0iCiAgIHZpZXdCb3g9IjAgMCA1Mi40NjY4NjIgNjEuNjM2MTkyIgogICB2ZXJzaW9uPSIxLjEiCiAgIGlkPSJzdmc1IgogICBpbmtzY2FwZTp2ZXJzaW9uPSIxLjEgKGM2OGUyMmMzODcsIDIwMjEtMDUtMjMpIgogICBzb2RpcG9kaTpkb2NuYW1lPSJDb25kZW5zb3Iuc3ZnIgogICB4bWxuczppbmtzY2FwZT0iaHR0cDovL3d3dy5pbmtzY2FwZS5vcmcvbmFtZXNwYWNlcy9pbmtzY2FwZSIKICAgeG1sbnM6c29kaXBvZGk9Imh0dHA6Ly9zb2RpcG9kaS5zb3VyY2Vmb3JnZS5uZXQvRFREL3NvZGlwb2RpLTAuZHRkIgogICB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciCiAgIHhtbG5zOnN2Zz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciPgogIDxzb2RpcG9kaTpuYW1lZHZpZXcKICAgICBpZD0ibmFtZWR2aWV3NyIKICAgICBwYWdlY29sb3I9IiNmZmZmZmYiCiAgICAgYm9yZGVyY29sb3I9IiM2NjY2NjYiCiAgICAgYm9yZGVyb3BhY2l0eT0iMS4wIgogICAgIGlua3NjYXBlOnBhZ2VzaGFkb3c9IjIiCiAgICAgaW5rc2NhcGU6cGFnZW9wYWNpdHk9IjAuMCIKICAgICBpbmtzY2FwZTpwYWdlY2hlY2tlcmJvYXJkPSIwIgogICAgIGlua3NjYXBlOmRvY3VtZW50LXVuaXRzPSJtbSIKICAgICBzaG93Z3JpZD0iZmFsc2UiCiAgICAgaW5rc2NhcGU6em9vbT0iMC44MDA4NzY4MiIKICAgICBpbmtzY2FwZTpjeD0iMTk0LjE2MjE5IgogICAgIGlua3NjYXBlOmN5PSIzMTAuMjg0OTIiCiAgICAgaW5rc2NhcGU6d2luZG93LXdpZHRoPSIyNTYwIgogICAgIGlua3NjYXBlOndpbmRvdy1oZWlnaHQ9IjEzNzciCiAgICAgaW5rc2NhcGU6d2luZG93LXg9Ii04IgogICAgIGlua3NjYXBlOndpbmRvdy15PSItOCIKICAgICBpbmtzY2FwZTp3aW5kb3ctbWF4aW1pemVkPSIxIgogICAgIGlua3NjYXBlOmN1cnJlbnQtbGF5ZXI9ImxheWVyMSIgLz4KICA8ZGVmcwogICAgIGlkPSJkZWZzMiIgLz4KICA8ZwogICAgIGlua3NjYXBlOmxhYmVsPSJMYXllciAxIgogICAgIGlua3NjYXBlOmdyb3VwbW9kZT0ibGF5ZXIiCiAgICAgaWQ9ImxheWVyMSIKICAgICB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMC40MjEwNzY0NSwtMC4wNTQzMzQzOSkiPgogICAgPGcKICAgICAgIHN0eWxlPSJmaWxsOm5vbmUiCiAgICAgICBpZD0iZzEwMjMiCiAgICAgICB0cmFuc2Zvcm09Im1hdHJpeCgwLjkwNDYwMTA3LDAsMCwwLjkwNDYwMTA3LC0yLjI5MjcyNjcsLTguNzQ3NDE1OSkiPgogICAgICA8cGF0aAogICAgICAgICBkPSJtIDU3LjQ2LDI4LjM1IGMgLTAuNTIwNCwtMC4wODkzIC0xLjA1MzksLTAuMDY1MyAtMS41NjQyLDAuMDcwNCAtMC41MTAyLDAuMTM1NyAtMC45ODUyLDAuMzc5OSAtMS4zOTI0LDAuNzE2IC0wLjQwNzMsMC4zMzYgLTAuNzM3MiwwLjc1NiAtMC45NjczLDEuMjMxMiAtMC4yMzAxLDAuNDc1MiAtMC4zNTUsMC45OTQ1IC0wLjM2NjEsMS41MjI0IHYgNC41NCBoIC0yLjU1IHYgLTE0IGggMy41NSBjIDAuMjY1MiwwIDAuNTE5NiwtMC4xMDU0IDAuNzA3MSwtMC4yOTI5IEMgNTUuMDY0NiwyMS45NDk2IDU1LjE3LDIxLjY5NTIgNTUuMTcsMjEuNDMgViAxNi4zIGMgMCwtMC4yNjUyIC0wLjEwNTQsLTAuNTE5NiAtMC4yOTI5LC0wLjcwNzEgQyA1NC42ODk2LDE1LjQwNTMgNTQuNDM1MiwxNS4zIDU0LjE3LDE1LjMgaCAtMy41NSB2IC00LjU3IGMgMCwtMC4yNjUyIC0wLjEwNTQsLTAuNTE5NiAtMC4yOTI5LC0wLjcwNzEgQyA1MC4xMzk2LDkuODM1MzQgNDkuODg1Miw5LjcyOTk4IDQ5LjYyLDkuNzI5OTggSCA0NyBjIC0wLjI2NTIsMCAtMC41MTk2LDAuMTA1MzYgLTAuNzA3MSwwLjI5MjkyIEMgNDYuMTA1NCwxMC4yMTA0IDQ2LDEwLjQ2NDggNDYsMTAuNzMgdiAyLjM1IEggMTIuMTcgdiAtMi4zNSBjIDAsLTAuMjY1MiAtMC4xMDU0LC0wLjUxOTYgLTAuMjkyOSwtMC43MDcxIEMgMTEuNjg5Niw5LjgzNTM0IDExLjQzNTIsOS43Mjk5OCAxMS4xNyw5LjcyOTk4IEggOC41NiBjIC0wLjI2NTIyLDAgLTAuNTE5NTcsMC4xMDUzNiAtMC43MDcxMSwwLjI5MjkyIEMgNy42NjUzNiwxMC4yMTA0IDcuNTYsMTAuNDY0OCA3LjU2LDEwLjczIFYgMTUuMyBIIDQgQyAzLjczNDc4LDE1LjMgMy40ODA0MywxNS40MDUzIDMuMjkyODksMTUuNTkyOSAzLjEwNTM2LDE1Ljc4MDQgMywxNi4wMzQ4IDMsMTYuMyB2IDUuMTEgYyAwLDAuMjY1MiAwLjEwNTM2LDAuNTE5NiAwLjI5Mjg5LDAuNzA3MSBDIDMuNDgwNDMsMjIuMzA0NiAzLjczNDc4LDIyLjQxIDQsMjIuNDEgSCA3LjU2IFYgNDEuNTkgSCA0IGMgLTAuMjY1MjIsMCAtMC41MTk1NywwLjEwNTMgLTAuNzA3MTEsMC4yOTI5IEMgMy4xMDUzNiw0Mi4wNzA0IDMsNDIuMzI0OCAzLDQyLjU5IHYgNS4xMSBjIDAsMC4yNjUyIDAuMTA1MzYsMC41MTk2IDAuMjkyODksMC43MDcxIEMgMy40ODA0Myw0OC41OTQ2IDMuNzM0NzgsNDguNyA0LDQ4LjcgaCAzLjU2IHYgNC41NyBjIDAsMC4yNjUyIDAuMTA1MzYsMC41MTk2IDAuMjkyODksMC43MDcxIEMgOC4wNDA0Myw1NC4xNjQ2IDguMjk0NzgsNTQuMjcgOC41Niw1NC4yNyBoIDIuNjEgYyAwLjI2NTIsMCAwLjUxOTYsLTAuMTA1NCAwLjcwNzEsLTAuMjkyOSBDIDEyLjA2NDYsNTMuNzg5NiAxMi4xNyw1My41MzUyIDEyLjE3LDUzLjI3IFYgNTAuOTIgSCA0NiB2IDIuMzUgYyAwLDAuMjY1MiAwLjEwNTQsMC41MTk2IDAuMjkyOSwwLjcwNzEgQyA0Ni40ODA0LDU0LjE2NDYgNDYuNzM0OCw1NC4yNyA0Nyw1NC4yNyBoIDIuNjIgYyAwLjI2NTIsMCAwLjUxOTYsLTAuMTA1NCAwLjcwNzEsLTAuMjkyOSBDIDUwLjUxNDYsNTMuNzg5NiA1MC42Miw1My41MzUyIDUwLjYyLDUzLjI3IFYgNDguNyBoIDMuNTUgYyAwLjI2NTIsMCAwLjUxOTYsLTAuMTA1NCAwLjcwNzEsLTAuMjkyOSBDIDU1LjA2NDYsNDguMjE5NiA1NS4xNyw0Ny45NjUyIDU1LjE3LDQ3LjcgdiAtMS41NiBjIDAuNzE5NiwwLjMxNTIgMS41MDg0LDAuNDM5MiAyLjI5LDAuMzYgMC45MzIsMCAxLjgyNjQsLTAuMzY3NCAyLjQ4OTEsLTEuMDIyNyBDIDYwLjYxMTksNDQuODIyMSA2MC45ODk1LDQzLjkzMTkgNjEsNDMgViAzMS44OSBDIDYwLjk5NzQsMzAuOTUxOSA2MC42MjM2LDMwLjA1MyA1OS45NjAyLDI5LjM4OTcgNTkuMjk2OSwyOC43MjY0IDU4LjM5ODEsMjguMzUyNiA1Ny40NiwyOC4zNSBaIG0gLTQuMjksMTAuMDggdiAzLjE2IGggLTIuNTUgdiAtMy4xNiB6IG0gMCwtMjEuMTMgdiAzLjExIEggNTAuNjIgViAxNy4zIFogbSAtMzIsLTIuMjIgaCAyLjQ4IHYgMzMuODQgaCAtMi41MiB6IG0gLTIsMzMuODQgSCAxNi42NSBWIDE1LjA4IGggMi40OCB6IG0gNi40OCwtMzMuODQgaCAyLjQ3IHYgMzMuODQgaCAtMi41MSB6IG0gNC40NywwIGggMi40OCB2IDMzLjg0IGggLTIuNTIgeiBtIDQuNDgsMCBIIDM3IHYgMzMuODQgaCAtMi40NCB6IG0gNC40OCwwIGggMi40OCBWIDQ4LjkyIEggMzkgWiBNIDUsMjAuNDEgViAxNy4zIGggMi41NiB2IDMuMTEgeiBNIDUsNDYuNyB2IC0zLjExIGggMi41NiB2IDMuMTEgeiBtIDUuMTcsNS41NyBIIDkuNTYgViAxMS43MyBoIDAuNjEgeiBtIDIsLTM3LjE5IGggMi40OCBWIDQ4LjkyIEggMTIuMTcgWiBNIDQzLjUyLDQ4LjkyIFYgMTUuMDggSCA0NiB2IDMzLjg0IHogbSA1LjEsMy4zNSBIIDQ4IFYgMTEuNzMgaCAwLjYyIHogbSA0LjU1LC01LjU3IGggLTIuNTUgdiAtMy4xMSBoIDIuNTUgeiBNIDU5LDQzIGMgLTAuMDAyNiwwLjQwNjcgLTAuMTY2MSwwLjc5NTggLTAuNDU0NiwxLjA4MjUgQyA1OC4yNTY5LDQ0LjM2OTEgNTcuODY2Nyw0NC41MyA1Ny40Niw0NC41MyBoIC0wLjc1IGMgLTAuNDA2NywwIC0wLjc5NjksLTAuMTYwOSAtMS4wODU0LC0wLjQ0NzUgQyA1NS4zMzYxLDQzLjc5NTggNTUuMTcyNiw0My40MDY3IDU1LjE3LDQzIFYgMzEuODkgYyAwLC0wLjQwODUgMC4xNjIyLC0wLjgwMDIgMC40NTExLC0xLjA4OSAwLjI4ODgsLTAuMjg4OCAwLjY4MDUsLTAuNDUxIDEuMDg4OSwtMC40NTEgaCAwLjc1IGMgMC40MDg0LDAgMC44MDAxLDAuMTYyMiAxLjA4ODksMC40NTEgQyA1OC44Mzc3LDMxLjA4OTggNTksMzEuNDgxNSA1OSwzMS44OSBaIgogICAgICAgICBmaWxsPSIjMDAwMDAwIgogICAgICAgICBpZD0icGF0aDEwMTQiIC8+CiAgICA8L2c+CiAgICA8ZwogICAgICAgaWQ9ImcxNDgyIgogICAgICAgdHJhbnNmb3JtPSJtYXRyaXgoMC4wMzY1NjEzNSwwLDAsMC4wMzY1NjEzNSwxNi41MDIwNjEsMzkuOTY5ODY2KSI+CiAgICAgIDxwYXRoCiAgICAgICAgIGQ9Ik0gNDQ3LjAyMzQ0LDg2LjY3OTY4OCA0NzcuNjk1MzEsNTYuMDAzOTA2IDM4OC42NDA2MiwzOC4wNzgxMjUgQyAzNDcuODM5ODQsMTMuMTY0MDYyIDMwMi4xNzk2OSwwIDI1Ni40MDIzNCwwIDE5NC4xMjUsMCAxMzQuMDE5NTMsMjMuNTc0MjE5IDg3LjE3NTc4MSw2NS40NzY1NjIgTCA1Ni41MDM5MDYsMzQuODAwNzgxIDM4LjU3NDIxOSwxMjMuODU5MzggQyAxMy42NjAxNTYsMTY0LjY2MDE2IDAuNSwyMTAuMzIwMzEgMC41LDI1Ni4wOTc2NiBjIDAsNjIuMjczNDMgMjMuNTc0MjE5LDEyMi4zNzg5IDY1LjQ3MjY1NiwxNjkuMjIyNjUgbCAtMzAuNjcxODc1LDMwLjY3MTg4IDg5LjA1ODU4OSwxNy45Mjk2OSBjIDQwLjgwMDc5LDI0LjkxNDA2IDg2LjQ2MDk0LDM4LjA3NDIxIDEzMi4yMzgyOSwzOC4wNzQyMSA2Mi4yNzM0MywwIDEyMi4zNzg5LC0yMy41NzQyMSAxNjkuMjIyNjUsLTY1LjQ3MjY1IGwgMzAuNzM0MzgsMzAuNzM0MzcgMTcuMzAwNzgsLTg3LjA4OTg0IDAuNDQ5MjIsLTAuNDQ5MjIgMC4zNTkzNywtMS44NDM3NSBjIDI0LjczODI4LC00MC41OTc2NiAzNy44MTI1LC04Ni4wOTM3NSAzNy44MTI1LC0xMzEuNzY1NjIgMCwtNjMuNTY2NDEgLTI0LjQxMDE1LC0xMjMuNTIzNDQgLTY1LjQ1MzEyLC0xNjkuNDI5NjkyIHogTSA3NC44NTkzNzUsOTUuNTYyNSA4Ni44ODI4MTIsMTA3LjU4NTk0IDk3LjQ4NDM3NSw5Ni45ODQzNzUgYyA0My4yMDcwMzUsLTQzLjIwNzAzMSA5OS42NDQ1MzUsLTY3IDE1OC45MTc5NjUsLTY3IDIzLjMwNDY5LDAgNDYuNjA5MzgsMy45MDIzNDQgNjkuMTA5MzgsMTEuNDY0ODQ0IGwgNi45Njg3NSwzNC44Mzk4NDMgQyAzMDguMTE3MTksNjUuNjc1NzgxIDI4Mi4zMTY0MSw2MC4xNTIzNDQgMjU2LjQ3NjU2LDYwLjE1MjM0NCBjIC00OS44Mzk4NCwwIC05OS43MDMxMiwxOS45NDUzMTIgLTEzNy43ODkwNiw1OC4wMzUxNTYgbCAtMTAuNjAxNTYsMTAuNTk3NjYgMTEuOTkyMTgsMTEuOTk2MDkgLTU2LjUzNTE1MSwxMS4zMDg1OSBDIDc3Ljg4NjcxOSw4MC40NDE0MDYgNzQuMjYxNzE5LDk4LjUzOTA2MiA3NC44NTkzNzUsOTUuNTYyNSBaIE0gOTYuMDU4NTk0LDQzNy42NDA2MiAxMDguMDg1OTQsNDI1LjYxNzE5IDk3LjQ4NDM3NSw0MTUuMDE1NjIgYyAtNDMuMjA3MDMxLC00My4yMDcwMyAtNjcsLTk5LjY0NDUzIC02NywtMTU4LjkxNzk2IDAsLTIzLjMwNDY5IDMuOTAyMzQ0LC00Ni42MDkzOCAxMS40NjA5MzcsLTY5LjEwOTM4IGwgMzQuODQzNzUsLTYuOTY4NzUgYyAtMTAuNjEzMjgxLDI0LjM2MzI4IC0xNi4xMzY3MTgsNTAuMTY0MDYgLTE2LjEzNjcxOCw3Ni4wMDM5MSAwLDQ5LjgzOTg0IDE5Ljk0NTMxMiw5OS43MDMxMiA1OC4wMzEyNDYsMTM3Ljc4OTA2IGwgMTAuNjAxNTcsMTAuNjAxNTYgMTEuOTk2MDksLTExLjk5MjE4IDExLjMwNDY5LDU2LjUzNTE1IEMgODAuOTM3NSw0MzQuNjEzMjggOTkuMDM5MDYyLDQzOC4yMzgyOCA5Ni4wNTg1OTQsNDM3LjY0MDYyIFogTSAxNTkuNjk1MzEsMzMxLjYwMTU2IDEyOS43MzQzNywzNjEuNTYyNSBDIDEwNC40MjE4OCwzMzEuNTI3MzQgOTAuNjM2NzE5LDI5NC41NzAzMSA5MC42MzY3MTksMjU2LjAxOTUzIGMgMCwtMjguNjY0MDYgOC4xNDQ1MzEsLTU3LjQwMjM0IDIzLjYwOTM3MSwtODMuNDkyMTkgbCA2Ni42NTIzNSwtMTMuMzI4MTIgLTI5Ljk2MDk0LC0yOS45NjA5NCBjIDMwLjAzMTI1LC0yNS4zMTY0IDY2Ljk4ODI4LC0zOS4xMDE1NjEgMTA1LjUzOTA2LC0zOS4xMDE1NjEgMjguNjY0MDYsMCA1Ny40MDIzNSw4LjE0ODQzNyA4My40OTIxOSwyMy42MDkzNzEgbCAxMy4zMzIwMyw2Ni42NTIzNSAyOS45NjA5NCwtMjkuOTYwOTQgYyAyNS4zMTY0LDMwLjAzNTE2IDM5LjA5NzY1LDY2Ljk5MjE5IDM5LjA5NzY1LDEwNS41NDI5NyAwLDI4LjY2NDA2IC04LjE0NDUzLDU3LjQwMjM0IC0yMy42MDU0Niw4My40OTIxOSBsIC02Ni42NTYyNSwxMy4zMjgxMiAyOS45NjA5MywyOS45NjA5NCBjIC0zMC4wMzEyNSwyNS4zMTY0IC02Ni45ODgyOCwzOS4xMDE1NiAtMTA1LjUzOTA2LDM5LjEwMTU2IC0yOC42NjQwNiwwIC01Ny40MDIzNCwtOC4xNDg0NCAtODMuNDkyMTksLTIzLjYwOTM3IHogbSAyNzguMzgyODEsODQuNzc3MzUgLTExLjk2NDg0LC0xMS45NjQ4NSAtMTAuNTk3NjYsMTAuNjAxNTcgYyAtNDMuMjA3MDMsNDMuMjAzMTIgLTk5LjY0NDUzLDY3IC0xNTguOTE3OTYsNjcgLTIzLjMwNDY5LDAgLTQ2LjYxMzI4LC0zLjkwNjI1IC02OS4xMTMyOCwtMTEuNDY0ODUgbCAtNi45NjQ4NSwtMzQuODM5ODQgYyAyNC4zNTkzOCwxMC42MTMyOCA1MC4xNjQwNiwxNi4xMzY3MiA3NiwxNi4xMzY3MiA0OS44NDM3NSwwIDk5LjcwMzEzLC0xOS45NDUzMiAxMzcuNzkyOTcsLTU4LjAzNTE2IGwgMTAuNjAxNTYsLTEwLjYwMTU2IC0xMS45OTYwOSwtMTEuOTkyMTkgNTYuNjI4OSwtMTEuMzI4MTMgYyAtNC4xNTYyNSwyMC41MjczNSAzLjEyODkxLC0xNC4yMzgyOCAtMTEuNDY4NzUsNTYuNDg4MjkgeiBNIDQ3MS4xMTMyOCwzMjUgbCAtMzQuOTAyMzQsNi45ODA0NyBjIDEwLjYwOTM3LC0yNC4zNjMyOCAxNi4xMzY3MiwtNTAuMTY0MDYgMTYuMTM2NzIsLTc2LjAwMzkxIDAsLTQ5LjczNDM3IC0xOS44NzExLC05OS42MjUgLTU4LjAzNTE2LC0xMzcuNzg5MDYgbCAtMTAuNjAxNTYsLTEwLjYwMTU2IC0xMS45OTIxOSwxMS45OTIxOSAtMTEuMzA4NTksLTU2LjUzNTE2MSBjIDcxLjY0ODQzLDE0LjMzOTg0MyA1My41NTA3OCwxMC43MTg3NSA1Ni41MzEyNSwxMS4zMTY0MDYgbCAtMTIuMDI3MzUsMTIuMDIzNDM3IDEwLjYwMTU3LDEwLjYwMTU2MyBjIDQzLjE5MTQsNDMuMTkxNDA1IDY2Ljk3NjU2LDk5LjcwMzEyNSA2Ni45NzY1NiwxNTkuMTI4OTA1IDAsMjMuMjQ2MDkgLTMuODc1LDQ2LjQ4MDQ3IC0xMS4zNzg5MSw2OC44ODY3MiB6IG0gMCwwIgogICAgICAgICBpZD0icGF0aDE0NzAiIC8+CiAgICAgIDxwYXRoCiAgICAgICAgIGQ9Im0gMzQ2LjQ0OTIyLDI3MC45OTIxOSBjIDAsLTQ2Ljk1MzEzIC03NS4wNTA3OCwtMTQxLjMwMDc4IC03OC4yNDIxOSwtMTQ1LjI5Mjk3IEwgMjU2LjUsMTExLjA2MjUgMjQ0Ljc5Mjk3LDEyNS42OTkyMiBjIC0zLjE5NTMxLDMuOTkyMTkgLTc4LjI0NjEsOTguMzM5ODQgLTc4LjI0NjEsMTQ1LjI5Mjk3IDAsNDkuNTk3NjUgNDAuMzUxNTcsODkuOTUzMTIgODkuOTUzMTMsODkuOTUzMTIgNDkuNTk3NjYsMCA4OS45NDkyMiwtNDAuMzU1NDcgODkuOTQ5MjIsLTg5Ljk1MzEyIHogTSAyNTYuNSwzMzAuOTYwOTQgYyAtMzMuMDY2NDEsMCAtNTkuOTY4NzUsLTI2LjkwMjM1IC01OS45Njg3NSwtNTkuOTY4NzUgMCwtMjEuOTgwNDcgMzIuNDg0MzcsLTc0LjMzOTg1IDU5Ljk2ODc1LC0xMTEuMzQzNzUgMjcuNDgwNDcsMzcuMDAzOSA1OS45NjQ4NCw4OS4zNjMyOCA1OS45NjQ4NCwxMTEuMzQzNzUgMCwzMy4wNjY0IC0yNi44OTg0Myw1OS45Njg3NSAtNTkuOTY0ODQsNTkuOTY4NzUgeiBtIDAsMCIKICAgICAgICAgaWQ9InBhdGgxNDcyIiAvPgogICAgPC9nPgogIDwvZz4KPC9zdmc+Cg==",
              fileName="//nas.ads.mwn.de/ge43kok/Desktop/SVGs/Condensor.svg"),
              Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              lineThickness=1)}),                                    Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Condensor_R134a;

    model Evaporator_R134a

      Modelica.Blocks.Interfaces.RealInput H_refri_in               "Inlet refrigerant specific enthalpy"    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-110,20})));
      Modelica.Blocks.Interfaces.RealInput MassFlow_refri_in        "Inlet refrigerant massflow"    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-110,-20})));
      Modelica.Blocks.Interfaces.RealInput MassFlow_refri_out       "Outlet refrigerant massflow"    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-110,-60})));
      Modelica.Blocks.Interfaces.RealOutput H_refri_out             "Outlet refrigerant specific enthalpy"    annotation (Placement(transformation(extent={{100,30},{120,50}})));
      Modelica.Blocks.Interfaces.RealInput T_air                    "Inlet air temperature" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-110,60})));
      Modelica.Blocks.Interfaces.RealOutput Pe                      "Refrigerant pressure in the condensor"   annotation (Placement(transformation(extent={{100,70},{120,90}})));
      Modelica.Blocks.Interfaces.RealOutput D_refri_out             "Outlet refrigerant density"   annotation (Placement(transformation(extent={{100,-10},{120,10}})));
      Modelica.Blocks.Interfaces.RealOutput Power_heat              "Heat power"   annotation (Placement(transformation(extent={{10,-10},{-10,10}},
            rotation=90,
            origin={0,-110})));
      Modelica.Blocks.Interfaces.RealOutput Superheat               "Temperautre level above the saturation threshold"   annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=90,
            origin={0,110})));
      Modelica.Blocks.Interfaces.RealOutput T_refri_out             "Outlet refrigerant temperauture" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={110,-40})));

    import si = Modelica.SIunits;

    /************************* Define all parameters ***************************************/
            parameter Real n = 1                                       "Number of tubes"                                                  annotation (Dialog(tab="Physical parameters"));
            parameter si.Length Ro = 9.5e-3                            "Outer radius of the tube"                                         annotation (Dialog(tab="Physical parameters"));
            parameter si.Length Ri = 8.25e-3                           "Inner radius of the tube"                                         annotation (Dialog(tab="Physical parameters"));
            parameter si.Length Lt = 15                                "Total length of each tube"                                        annotation (Dialog(tab="Physical parameters"));
            parameter si.CoefficientOfHeatTransfer Alpa_i_vapor = 400  "Inner heat transfer coefficient of the vapor zone"                annotation (Dialog(tab="Physical parameters"));
            parameter si.CoefficientOfHeatTransfer Alpa_i_satu = 2000  "Inner heat transfer coefficient of the two phase zone"            annotation (Dialog(tab="Physical parameters"));
            parameter si.CoefficientOfHeatTransfer Alpa_o = 300        "Outer heat transfer coefficient of evaporator"                    annotation (Dialog(tab="Physical parameters"));
            parameter si.SpecificHeatCapacity Cp_wall = 385            "Specific heat of pipe wall material"                              annotation (Dialog(tab="Physical parameters"));
            parameter si.Density D_wall = 8960                         "Density of pipe wall material"                                    annotation (Dialog(tab="Physical parameters"));
            parameter si.Pressure Pe_initial = 4e5                     "Initial refrigerant pressure in the evaporator"                   annotation (Dialog(tab="Initial conditions"));
            parameter si.SpecificEnthalpy H_refri_out_initial = 410e3  "Initial output refrigerant specific enthalpy of the evaporator"   annotation (Dialog(tab="Initial conditions"));
            parameter si.Temperature T_wall_1_initial = 293.15         "For evaporator, initial wall temperature of the two-phase zone"   annotation (Dialog(tab="Initial conditions"));
            parameter si.Temperature T_wall_2_initial = 293.15         "For evaporator, initial wall temperature of the vapor zone"       annotation (Dialog(tab="Initial conditions"));
            parameter si.Length L1_initial = 5                         "For evaporator, initial length of the two-phase zone"             annotation (Dialog(tab="Initial conditions"));
    /************************* Define all variables ***************************************/
            Real Gama                         "Average void fraction";
            Real a;

            si.Length L1                      "For evaporator, length of the two-phase zone";
            si.Temperature T_wall_1           "For evaporator, wall temperature of the two-phase zone";
            si.Temperature T_wall_2           "For evaporator, wall temperature of the vapor zone";
            si.Temperature T_wall_boundary    "Temperature at the boundary between two zones";
            si.MassFlowRate MassFlow_boundary "Massflow at the boundary between two zones";

            Modelica.Media.R134a.R134a_ph.ThermodynamicState Vapor_state  "Thermal dynamic state variable for refrigerant vapor";
            Modelica.Media.R134a.R134a_ph.SaturationProperties Satu_state "Thermal dynamic state variable for two-phase refrigerant mixture";

            si.Area A_cs_o "Tube outer crossection area (pipe thickness counted)";
            si.Area A_cs_i "Tube inner crossection area (hydraulic)";
            si.Area A_sf_i "Tube inner surface area (contact area between refrigerant and pipe wall)";
            si.Area A_sf_o "Tube outer surface area (contact area between air and pipe wall)";

            si.Temperature T_refri_satu             "Refrigerant temperature in two-phase zone";
            si.SpecificEnthalpy H_refri_satu_liquid "Refrigerant specific enthalpy for saturated liquid state";
            si.SpecificEnthalpy H_refri_satu_gas    "Refrigerant specific enthalpy for saturated gas state";
            si.Density D_refri_satu_liquid          "Refrigerant density for saturated liquid state";
            si.Density D_refri_satu_gas             "Refrigerant density for saturated gas state";

            si.SpecificEnthalpy H_refri_vapor_avg           "Refrigerant specific enthalpy for vapor zone (averaged)";
            si.Density D_refri_vapor_avg           "Refrigerant density for vapor zone (averaged)";
            si.Temperature T_refri_vapor_avg "Refrigerant temperature for vapor zone (averaged)";

            Real dDdH_refri_at_Pe_vapor                   "Refrigerant thermal dynamics property, calculated from data";
            Real dDdP_refri_at_H_vapor                    "Refrigerant thermal dynamics property, calculated from data";
            Real dDdP_refri_satu_liquid                   "Refrigerant thermal dynamics property, calculated from data";
            Real dDdP_refri_satu_gas                      "Refrigerant thermal dynamics property, calculated from data";
            Real dHdP_refri_satu_liquid                   "Refrigerant thermal dynamics property, calculated from data";
            Real dHdP_refri_satu_gas                      "Refrigerant thermal dynamics property, calculated from data";

    initial equation

            Pe = Pe_initial;
            H_refri_out = H_refri_out_initial;
            T_wall_1 = T_wall_1_initial;
            T_wall_2 = T_wall_2_initial;
            L1 = L1_initial;

    equation
            A_cs_o = Ro^2 * Modelica.Constants.pi "Outer cross-section area of condensor tubes";
            A_cs_i = Ri^2 * Modelica.Constants.pi "Inner cross-section area of condensor tubes";
            A_sf_i = n * Ri * 2 * Modelica.Constants.pi * Lt "Inner heat transfer effective area ";
            A_sf_o = n * Ro * 2 * Modelica.Constants.pi * Lt "Outer heat transfer effective area ";

            a = (D_refri_satu_gas/D_refri_satu_liquid)^(2/3);
            Gama = (( 1 - a * (1 - log(a)))   / (1-a)^2);

    /************************* Fluid properties calculation for the two-phase zone ***************/
            T_refri_satu = Modelica.Media.R134a.R134a_ph.saturationTemperature(Pe);
            Satu_state.psat = Pe;
            Satu_state.Tsat = T_refri_satu;
            dDdP_refri_satu_liquid = Modelica.Media.R134a.R134a_ph.dBubbleDensity_dPressure(Satu_state);
            dDdP_refri_satu_gas = Modelica.Media.R134a.R134a_ph.dDewDensity_dPressure(Satu_state);
            dHdP_refri_satu_liquid  = Modelica.Media.R134a.R134a_ph.dBubbleEnthalpy_dPressure(Satu_state);
            dHdP_refri_satu_gas  = Modelica.Media.R134a.R134a_ph.dDewEnthalpy_dPressure(Satu_state);
            H_refri_satu_liquid = Modelica.Media.R134a.R134a_ph.bubbleEnthalpy(Satu_state);
            H_refri_satu_gas = Modelica.Media.R134a.R134a_ph.dewEnthalpy(Satu_state);
            D_refri_satu_liquid = Modelica.Media.R134a.R134a_ph.bubbleDensity(Satu_state);
            D_refri_satu_gas = Modelica.Media.R134a.R134a_ph.dewDensity(Satu_state);

    /************************* Fluid properties calculation for the vapor zone ***************/
            Vapor_state.phase = 1;
            Vapor_state.h = H_refri_vapor_avg;
            Vapor_state.d = D_refri_vapor_avg;
            Vapor_state.T = T_refri_vapor_avg;
            Vapor_state.p = Pe;
            D_refri_vapor_avg = Modelica.Media.R134a.R134a_ph.density_ph(Pe, H_refri_vapor_avg, 1);
            T_refri_vapor_avg = Modelica.Media.R134a.R134a_ph.temperature_ph(Pe, H_refri_vapor_avg, 1);
            dDdH_refri_at_Pe_vapor = Modelica.Media.R134a.R134a_ph.density_derh_p(Vapor_state);
            dDdP_refri_at_H_vapor =  Modelica.Media.R134a.R134a_ph.density_derp_h(Vapor_state);
            D_refri_out =Modelica.Media.R134a.R134a_ph.density_ph(Pe,H_refri_out,1);

            T_wall_boundary = T_wall_1                                   "The boundary temperature between two-phase and vapor zone is set to be the two-phase zone temperature";
            H_refri_vapor_avg = 0.5 * (H_refri_satu_gas + H_refri_out)   "The refrigerant specific enthalpy in vapor zone is an averaged valve between input and output to the zone";
            T_refri_out = 2 * T_refri_vapor_avg - T_refri_satu;
            Superheat = 2 * (T_refri_vapor_avg - T_refri_satu);
            Power_heat = Alpa_o * A_sf_o * (L1/Lt) * (T_air - T_wall_1) + Alpa_o * A_sf_o * ((Lt - L1)/Lt) * (T_air - T_wall_2);

    /************************* "Conservation of mass in the two-phase region" *****************************************/
            (dDdP_refri_satu_liquid * (1 - Gama) + dDdP_refri_satu_gas * Gama) * A_cs_i * L1 * der(Pe)
            + (D_refri_satu_liquid - D_refri_satu_gas) * (1 - Gama) * A_cs_i * der(L1)
            =
            MassFlow_refri_in - MassFlow_boundary;
            //+ (D_refri_satu_gas - D_refri_satu_liquid) * A_cs_i * L1 * der(Gama)

     /************************* "Conservation of mass in the vapor region" ********************************************/
            (dDdP_refri_at_H_vapor + 0.5 * dDdH_refri_at_Pe_vapor * dHdP_refri_satu_gas) * A_cs_i * (Lt - L1) * der(Pe)
            + 0.5*dDdH_refri_at_Pe_vapor*A_cs_i*(Lt - L1) * der(H_refri_out)
            + (D_refri_satu_gas - D_refri_vapor_avg) * A_cs_i * der(L1)
            =
            MassFlow_boundary - MassFlow_refri_out;

     /************************* "Conservation of energy in the two-phase region" **************************************/
            ((dDdP_refri_satu_liquid * H_refri_satu_liquid + dHdP_refri_satu_liquid  * D_refri_satu_liquid) * (1 - Gama) + (dDdP_refri_satu_gas * H_refri_satu_gas + dHdP_refri_satu_gas  * D_refri_satu_gas) * Gama - 1) * A_cs_i * L1 * der(Pe)
            + (1 - Gama) * (D_refri_satu_liquid * H_refri_satu_liquid - D_refri_satu_gas * H_refri_satu_gas) * A_cs_i * der(L1)
            =
            MassFlow_refri_in * H_refri_in
            - MassFlow_boundary * H_refri_satu_gas
            + Alpa_i_satu * A_sf_i * (L1/Lt) * (T_wall_1 - T_refri_satu);
            //+ (D_refri_satu_gas * H_refri_satu_gas - D_refri_satu_liquid * H_refri_satu_liquid) * A_cs_i * L1 * der(Gama)

    /************************* "Conservation of energy in the vapor region" ******************************************/
            ((dDdP_refri_at_H_vapor + 0.5 * dHdP_refri_satu_gas *dDdH_refri_at_Pe_vapor) * H_refri_vapor_avg + 0.5 * dHdP_refri_satu_gas *D_refri_vapor_avg - 1) * A_cs_i * (Lt - L1) * der(Pe)
            + 0.5 * (dDdH_refri_at_Pe_vapor*H_refri_vapor_avg + D_refri_vapor_avg) * A_cs_i * (Lt - L1) * der(H_refri_out)
            + (D_refri_satu_gas*H_refri_satu_gas - D_refri_vapor_avg*H_refri_vapor_avg) * A_cs_i*der(L1)
            =
            MassFlow_boundary * H_refri_satu_gas
            - MassFlow_refri_out * H_refri_out
            + Alpa_i_vapor * A_sf_i * ((Lt - L1)/Lt) * (T_wall_2 - T_refri_vapor_avg);

    /************************* "Conservation of energy of pipe wall (two-phase region)"********************************/
            Cp_wall * D_wall * (A_cs_o - A_cs_i) * L1 * der(T_wall_1)
            =
            Alpa_i_satu * A_sf_i * (L1/Lt) * (T_refri_satu - T_wall_1)
            + Alpa_o * A_sf_o * (L1/Lt) * (T_air - T_wall_1);

    /************************* "Conservation of energy of pipe wall (vapor region)"***********************************/
            Cp_wall * D_wall * (A_cs_o - A_cs_i) * (Lt - L1) * (der(T_wall_2) - ((T_wall_2 - T_wall_boundary)/(Lt - L1)) * der(L1))
            =
            Alpa_i_vapor * A_sf_i * ((Lt - L1)/Lt) * ( T_refri_vapor_avg - T_wall_2)
            + Alpa_o * A_sf_o * ((Lt - L1)/Lt) * (T_air - T_wall_2);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Bitmap(
              extent={{-98,-98},{98,98}},
              imageSource="PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiIHN0YW5kYWxvbmU9Im5vIj8+CjwhLS0gQ3JlYXRlZCB3aXRoIElua3NjYXBlIChodHRwOi8vd3d3Lmlua3NjYXBlLm9yZy8pIC0tPgoKPHN2ZwogICB3aWR0aD0iNTIuNDY2ODYybW0iCiAgIGhlaWdodD0iNjEuNjM2MTkybW0iCiAgIHZpZXdCb3g9IjAgMCA1Mi40NjY4NjIgNjEuNjM2MTkyIgogICB2ZXJzaW9uPSIxLjEiCiAgIGlkPSJzdmc1IgogICBpbmtzY2FwZTp2ZXJzaW9uPSIxLjEgKGM2OGUyMmMzODcsIDIwMjEtMDUtMjMpIgogICBzb2RpcG9kaTpkb2NuYW1lPSJFdmFwb3JhdG9yLnN2ZyIKICAgeG1sbnM6aW5rc2NhcGU9Imh0dHA6Ly93d3cuaW5rc2NhcGUub3JnL25hbWVzcGFjZXMvaW5rc2NhcGUiCiAgIHhtbG5zOnNvZGlwb2RpPSJodHRwOi8vc29kaXBvZGkuc291cmNlZm9yZ2UubmV0L0RURC9zb2RpcG9kaS0wLmR0ZCIKICAgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIgogICB4bWxuczpzdmc9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KICA8c29kaXBvZGk6bmFtZWR2aWV3CiAgICAgaWQ9Im5hbWVkdmlldzciCiAgICAgcGFnZWNvbG9yPSIjZmZmZmZmIgogICAgIGJvcmRlcmNvbG9yPSIjNjY2NjY2IgogICAgIGJvcmRlcm9wYWNpdHk9IjEuMCIKICAgICBpbmtzY2FwZTpwYWdlc2hhZG93PSIyIgogICAgIGlua3NjYXBlOnBhZ2VvcGFjaXR5PSIwLjAiCiAgICAgaW5rc2NhcGU6cGFnZWNoZWNrZXJib2FyZD0iMCIKICAgICBpbmtzY2FwZTpkb2N1bWVudC11bml0cz0ibW0iCiAgICAgc2hvd2dyaWQ9ImZhbHNlIgogICAgIGlua3NjYXBlOnpvb209IjAuODAwODc2ODIiCiAgICAgaW5rc2NhcGU6Y3g9IjE5NC4xNjIxOSIKICAgICBpbmtzY2FwZTpjeT0iMzEwLjI4NDkyIgogICAgIGlua3NjYXBlOndpbmRvdy13aWR0aD0iMjU2MCIKICAgICBpbmtzY2FwZTp3aW5kb3ctaGVpZ2h0PSIxMzc3IgogICAgIGlua3NjYXBlOndpbmRvdy14PSItOCIKICAgICBpbmtzY2FwZTp3aW5kb3cteT0iLTgiCiAgICAgaW5rc2NhcGU6d2luZG93LW1heGltaXplZD0iMSIKICAgICBpbmtzY2FwZTpjdXJyZW50LWxheWVyPSJsYXllcjEiIC8+CiAgPGRlZnMKICAgICBpZD0iZGVmczIiIC8+CiAgPGcKICAgICBpbmtzY2FwZTpsYWJlbD0iTGF5ZXIgMSIKICAgICBpbmtzY2FwZTpncm91cG1vZGU9ImxheWVyIgogICAgIGlkPSJsYXllcjEiCiAgICAgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTAuNDIxMDc2NDUsLTAuMDU0MzM0MzkpIj4KICAgIDxnCiAgICAgICBpZD0iZzEyMTkiCiAgICAgICB0cmFuc2Zvcm09Im1hdHJpeCgzLjQxODk2NDcsMCwwLDMuNDE4OTY0NywtMjkxLjg2ODI1LC05Mi43NTE1MjkpIj4KICAgICAgPGcKICAgICAgICAgc3R5bGU9ImZpbGw6bm9uZSIKICAgICAgICAgaWQ9ImcxMDIzIgogICAgICAgICB0cmFuc2Zvcm09Im1hdHJpeCgwLjI2NDU4MzMzLDAsMCwwLjI2NDU4MzMzLDg0LjY5Njg0NSwyNC41NzAwNDQpIj4KICAgICAgICA8cGF0aAogICAgICAgICAgIGQ9Im0gNTcuNDYsMjguMzUgYyAtMC41MjA0LC0wLjA4OTMgLTEuMDUzOSwtMC4wNjUzIC0xLjU2NDIsMC4wNzA0IC0wLjUxMDIsMC4xMzU3IC0wLjk4NTIsMC4zNzk5IC0xLjM5MjQsMC43MTYgLTAuNDA3MywwLjMzNiAtMC43MzcyLDAuNzU2IC0wLjk2NzMsMS4yMzEyIC0wLjIzMDEsMC40NzUyIC0wLjM1NSwwLjk5NDUgLTAuMzY2MSwxLjUyMjQgdiA0LjU0IGggLTIuNTUgdiAtMTQgaCAzLjU1IGMgMC4yNjUyLDAgMC41MTk2LC0wLjEwNTQgMC43MDcxLC0wLjI5MjkgQyA1NS4wNjQ2LDIxLjk0OTYgNTUuMTcsMjEuNjk1MiA1NS4xNywyMS40MyBWIDE2LjMgYyAwLC0wLjI2NTIgLTAuMTA1NCwtMC41MTk2IC0wLjI5MjksLTAuNzA3MSBDIDU0LjY4OTYsMTUuNDA1MyA1NC40MzUyLDE1LjMgNTQuMTcsMTUuMyBoIC0zLjU1IHYgLTQuNTcgYyAwLC0wLjI2NTIgLTAuMTA1NCwtMC41MTk2IC0wLjI5MjksLTAuNzA3MSBDIDUwLjEzOTYsOS44MzUzNCA0OS44ODUyLDkuNzI5OTggNDkuNjIsOS43Mjk5OCBIIDQ3IGMgLTAuMjY1MiwwIC0wLjUxOTYsMC4xMDUzNiAtMC43MDcxLDAuMjkyOTIgQyA0Ni4xMDU0LDEwLjIxMDQgNDYsMTAuNDY0OCA0NiwxMC43MyB2IDIuMzUgSCAxMi4xNyB2IC0yLjM1IGMgMCwtMC4yNjUyIC0wLjEwNTQsLTAuNTE5NiAtMC4yOTI5LC0wLjcwNzEgQyAxMS42ODk2LDkuODM1MzQgMTEuNDM1Miw5LjcyOTk4IDExLjE3LDkuNzI5OTggSCA4LjU2IGMgLTAuMjY1MjIsMCAtMC41MTk1NywwLjEwNTM2IC0wLjcwNzExLDAuMjkyOTIgQyA3LjY2NTM2LDEwLjIxMDQgNy41NiwxMC40NjQ4IDcuNTYsMTAuNzMgViAxNS4zIEggNCBDIDMuNzM0NzgsMTUuMyAzLjQ4MDQzLDE1LjQwNTMgMy4yOTI4OSwxNS41OTI5IDMuMTA1MzYsMTUuNzgwNCAzLDE2LjAzNDggMywxNi4zIHYgNS4xMSBjIDAsMC4yNjUyIDAuMTA1MzYsMC41MTk2IDAuMjkyODksMC43MDcxIEMgMy40ODA0MywyMi4zMDQ2IDMuNzM0NzgsMjIuNDEgNCwyMi40MSBIIDcuNTYgViA0MS41OSBIIDQgYyAtMC4yNjUyMiwwIC0wLjUxOTU3LDAuMTA1MyAtMC43MDcxMSwwLjI5MjkgQyAzLjEwNTM2LDQyLjA3MDQgMyw0Mi4zMjQ4IDMsNDIuNTkgdiA1LjExIGMgMCwwLjI2NTIgMC4xMDUzNiwwLjUxOTYgMC4yOTI4OSwwLjcwNzEgQyAzLjQ4MDQzLDQ4LjU5NDYgMy43MzQ3OCw0OC43IDQsNDguNyBoIDMuNTYgdiA0LjU3IGMgMCwwLjI2NTIgMC4xMDUzNiwwLjUxOTYgMC4yOTI4OSwwLjcwNzEgQyA4LjA0MDQzLDU0LjE2NDYgOC4yOTQ3OCw1NC4yNyA4LjU2LDU0LjI3IGggMi42MSBjIDAuMjY1MiwwIDAuNTE5NiwtMC4xMDU0IDAuNzA3MSwtMC4yOTI5IEMgMTIuMDY0Niw1My43ODk2IDEyLjE3LDUzLjUzNTIgMTIuMTcsNTMuMjcgViA1MC45MiBIIDQ2IHYgMi4zNSBjIDAsMC4yNjUyIDAuMTA1NCwwLjUxOTYgMC4yOTI5LDAuNzA3MSBDIDQ2LjQ4MDQsNTQuMTY0NiA0Ni43MzQ4LDU0LjI3IDQ3LDU0LjI3IGggMi42MiBjIDAuMjY1MiwwIDAuNTE5NiwtMC4xMDU0IDAuNzA3MSwtMC4yOTI5IEMgNTAuNTE0Niw1My43ODk2IDUwLjYyLDUzLjUzNTIgNTAuNjIsNTMuMjcgViA0OC43IGggMy41NSBjIDAuMjY1MiwwIDAuNTE5NiwtMC4xMDU0IDAuNzA3MSwtMC4yOTI5IEMgNTUuMDY0Niw0OC4yMTk2IDU1LjE3LDQ3Ljk2NTIgNTUuMTcsNDcuNyB2IC0xLjU2IGMgMC43MTk2LDAuMzE1MiAxLjUwODQsMC40MzkyIDIuMjksMC4zNiAwLjkzMiwwIDEuODI2NCwtMC4zNjc0IDIuNDg5MSwtMS4wMjI3IEMgNjAuNjExOSw0NC44MjIxIDYwLjk4OTUsNDMuOTMxOSA2MSw0MyBWIDMxLjg5IEMgNjAuOTk3NCwzMC45NTE5IDYwLjYyMzYsMzAuMDUzIDU5Ljk2MDIsMjkuMzg5NyA1OS4yOTY5LDI4LjcyNjQgNTguMzk4MSwyOC4zNTI2IDU3LjQ2LDI4LjM1IFogbSAtNC4yOSwxMC4wOCB2IDMuMTYgaCAtMi41NSB2IC0zLjE2IHogbSAwLC0yMS4xMyB2IDMuMTEgSCA1MC42MiBWIDE3LjMgWiBtIC0zMiwtMi4yMiBoIDIuNDggdiAzMy44NCBoIC0yLjUyIHogbSAtMiwzMy44NCBIIDE2LjY1IFYgMTUuMDggaCAyLjQ4IHogbSA2LjQ4LC0zMy44NCBoIDIuNDcgdiAzMy44NCBoIC0yLjUxIHogbSA0LjQ3LDAgaCAyLjQ4IHYgMzMuODQgaCAtMi41MiB6IG0gNC40OCwwIEggMzcgdiAzMy44NCBoIC0yLjQ0IHogbSA0LjQ4LDAgaCAyLjQ4IFYgNDguOTIgSCAzOSBaIE0gNSwyMC40MSBWIDE3LjMgaCAyLjU2IHYgMy4xMSB6IE0gNSw0Ni43IHYgLTMuMTEgaCAyLjU2IHYgMy4xMSB6IG0gNS4xNyw1LjU3IEggOS41NiBWIDExLjczIGggMC42MSB6IG0gMiwtMzcuMTkgaCAyLjQ4IFYgNDguOTIgSCAxMi4xNyBaIE0gNDMuNTIsNDguOTIgViAxNS4wOCBIIDQ2IHYgMzMuODQgeiBtIDUuMSwzLjM1IEggNDggViAxMS43MyBoIDAuNjIgeiBtIDQuNTUsLTUuNTcgaCAtMi41NSB2IC0zLjExIGggMi41NSB6IE0gNTksNDMgYyAtMC4wMDI2LDAuNDA2NyAtMC4xNjYxLDAuNzk1OCAtMC40NTQ2LDEuMDgyNSBDIDU4LjI1NjksNDQuMzY5MSA1Ny44NjY3LDQ0LjUzIDU3LjQ2LDQ0LjUzIGggLTAuNzUgYyAtMC40MDY3LDAgLTAuNzk2OSwtMC4xNjA5IC0xLjA4NTQsLTAuNDQ3NSBDIDU1LjMzNjEsNDMuNzk1OCA1NS4xNzI2LDQzLjQwNjcgNTUuMTcsNDMgViAzMS44OSBjIDAsLTAuNDA4NSAwLjE2MjIsLTAuODAwMiAwLjQ1MTEsLTEuMDg5IDAuMjg4OCwtMC4yODg4IDAuNjgwNSwtMC40NTEgMS4wODg5LC0wLjQ1MSBoIDAuNzUgYyAwLjQwODQsMCAwLjgwMDEsMC4xNjIyIDEuMDg4OSwwLjQ1MSBDIDU4LjgzNzcsMzEuMDg5OCA1OSwzMS40ODE1IDU5LDMxLjg5IFoiCiAgICAgICAgICAgZmlsbD0iIzAwMDAwMCIKICAgICAgICAgICBpZD0icGF0aDEwMTQiIC8+CiAgICAgIDwvZz4KICAgICAgPGcKICAgICAgICAgaWQ9ImcxMTU1IgogICAgICAgICB0cmFuc2Zvcm09Im1hdHJpeCgwLjAxNzk3NzQ0LDAsMCwwLjAxNzk3NzQ0LDg4Ljk5NzE0NSwzOS4xMzE3NSkiPgogICAgICAgIDxnCiAgICAgICAgICAgaWQ9ImcxMDgyIj4KCTxnCiAgIGlkPSJnMTA4MCI+CgkJPHBhdGgKICAgZD0ibSAyOTYsNDggYyAtMzkuNzA0LDAgLTcyLDMyLjMwNCAtNzIsNzIgMCw0LjQxNiAzLjU3Niw4IDgsOCA0LjQyNCwwIDgsLTMuNTg0IDgsLTggMCwtMzAuODggMjUuMTI4LC01NiA1NiwtNTYgMzAuODcyLDAgNTYsMjUuMTIgNTYsNTYgMCwzMC44OCAtMjUuMTI4LDU2IC01Niw1NiBIIDggYyAtNC40MTYsMCAtOCwzLjU4NCAtOCw4IDAsNC40MTYgMy41ODQsOCA4LDggaCAyODggYyAzOS43MDQsMCA3MiwtMzIuMzA0IDcyLC03MiAwLC0zOS42OTYgLTMyLjI5NiwtNzIgLTcyLC03MiB6IgogICBpZD0icGF0aDEwNzgiIC8+CgoJPC9nPgoKPC9nPgogICAgICAgIDxnCiAgICAgICAgICAgaWQ9ImcxMDg4Ij4KCTxnCiAgIGlkPSJnMTA4NiI+CgkJPHBhdGgKICAgZD0ibSAxNDQsMzIgYyAtMzAuODgsMCAtNTYsMjUuMTIgLTU2LDU2IDAsNC40MTYgMy41ODQsOCA4LDggNC40MTYsMCA4LC0zLjU4NCA4LC04IDAsLTIyLjA1NiAxNy45NDQsLTQwIDQwLC00MCAyMi4wNTYsMCA0MCwxNy45NDQgNDAsNDAgMCwyMi4wNTYgLTE3Ljk0NCw0MCAtNDAsNDAgSCA4IGMgLTQuNDE2LDAgLTgsMy41ODQgLTgsOCAwLDQuNDE2IDMuNTg0LDggOCw4IGggMTM2IGMgMzAuODgsMCA1NiwtMjUuMTIgNTYsLTU2IDAsLTMwLjg4IC0yNS4xMiwtNTYgLTU2LC01NiB6IgogICBpZD0icGF0aDEwODQiIC8+CgoJPC9nPgoKPC9nPgogICAgICAgIDxnCiAgICAgICAgICAgaWQ9ImcxMDk0Ij4KCTxnCiAgIGlkPSJnMTA5MiI+CgkJPHBhdGgKICAgZD0iTSAyODAsMjI0IEggOCBjIC00LjQxNiwwIC04LDMuNTg0IC04LDggMCw0LjQxNiAzLjU4NCw4IDgsOCBoIDI3MiBjIDIyLjA1NiwwIDQwLDE3Ljk0NCA0MCw0MCAwLDIyLjA1NiAtMTcuOTQ0LDQwIC00MCw0MCAtMjIuMDU2LDAgLTQwLC0xNy45NDQgLTQwLC00MCAwLC00LjQxNiAtMy41NzYsLTggLTgsLTggLTQuNDI0LDAgLTgsMy41ODQgLTgsOCAwLDMwLjg4IDI1LjEyOCw1NiA1Niw1NiAzMC44NzIsMCA1NiwtMjUuMTIgNTYsLTU2IDAsLTMwLjg4IC0yNS4xMjgsLTU2IC01NiwtNTYgeiIKICAgaWQ9InBhdGgxMDkwIiAvPgoKCTwvZz4KCjwvZz4KICAgICAgICA8ZwogICAgICAgICAgIGlkPSJnMTA5NiI+CjwvZz4KICAgICAgICA8ZwogICAgICAgICAgIGlkPSJnMTA5OCI+CjwvZz4KICAgICAgICA8ZwogICAgICAgICAgIGlkPSJnMTEwMCI+CjwvZz4KICAgICAgICA8ZwogICAgICAgICAgIGlkPSJnMTEwMiI+CjwvZz4KICAgICAgICA8ZwogICAgICAgICAgIGlkPSJnMTEwNCI+CjwvZz4KICAgICAgICA8ZwogICAgICAgICAgIGlkPSJnMTEwNiI+CjwvZz4KICAgICAgICA8ZwogICAgICAgICAgIGlkPSJnMTEwOCI+CjwvZz4KICAgICAgICA8ZwogICAgICAgICAgIGlkPSJnMTExMCI+CjwvZz4KICAgICAgICA8ZwogICAgICAgICAgIGlkPSJnMTExMiI+CjwvZz4KICAgICAgICA8ZwogICAgICAgICAgIGlkPSJnMTExNCI+CjwvZz4KICAgICAgICA8ZwogICAgICAgICAgIGlkPSJnMTExNiI+CjwvZz4KICAgICAgICA8ZwogICAgICAgICAgIGlkPSJnMTExOCI+CjwvZz4KICAgICAgICA8ZwogICAgICAgICAgIGlkPSJnMTEyMCI+CjwvZz4KICAgICAgICA8ZwogICAgICAgICAgIGlkPSJnMTEyMiI+CjwvZz4KICAgICAgICA8ZwogICAgICAgICAgIGlkPSJnMTEyNCI+CjwvZz4KICAgICAgPC9nPgogICAgPC9nPgogIDwvZz4KPC9zdmc+Cg==",
              fileName="//nas.ads.mwn.de/ge43kok/Desktop/SVGs/Evaporator.svg"),
              Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              lineThickness=1)}),                                    Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Evaporator_R134a;

    model Compressor

      parameter Modelica.SIunits.Volume V_c = 150e-6 "Effective displacement volume";
      parameter Real R = 81.49 "Gas constant of the refrigerant";
      parameter Real k = 1.2 "Empirical parameter related to the refrigerant";

      Real eta_v "Volumetric efficiency";

     Modelica.Blocks.Interfaces.RealInput Speed "Motor rotation speed, in rads/s"
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=270,
            origin={0,-110})));
      Modelica.Blocks.Interfaces.RealInput D_in "Density of the refrigerant at the inlet" annotation (Placement(
            transformation(extent={{-120,30},{-100,50}}), iconTransformation(
              extent={{-120,30},{-100,50}})));
      Modelica.Blocks.Interfaces.RealOutput MassFlow
        "Massflow output of the compressor"
        annotation (Placement(transformation(extent={{100,50},{120,70}}),
            iconTransformation(extent={{100,50},{120,70}})));
      Modelica.Blocks.Interfaces.RealInput H_in
        "Specific enthalpy at the compressor inlet" annotation (Placement(
            transformation(extent={{-120,-10},{-100,10}}), iconTransformation(
              extent={{-120,-10},{-100,10}})));
      Modelica.Blocks.Interfaces.RealOutput H_out
        "Specific enthalpy at the compressor outlet" annotation (Placement(
            transformation(extent={{100,10},{120,30}}), iconTransformation(extent=
               {{100,10},{120,30}})));
      Modelica.Blocks.Interfaces.RealInput Pc "Presure in the condensor"
        annotation (Placement(transformation(extent={{-120,-90},{-100,-70}}),
            iconTransformation(extent={{-120,-90},{-100,-70}})));
      Modelica.Blocks.Interfaces.RealInput Pe "Presure in the evaporator"
        annotation (Placement(transformation(extent={{-120,-50},{-100,-30}})));

      Modelica.Blocks.Interfaces.RealOutput Power_Shaft
        "Shaft power of the compressor" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={0,110}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={0,110})));
     Modelica.Blocks.Interfaces.RealInput T_in "Temperature of the refrigerant at the inlet"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-110,80})));
    equation

      eta_v = 1 - 0.043*Pc/Pe;

      MassFlow = (Speed/6.2832) * V_c *  eta_v * D_in;

      Power_Shaft = MassFlow*(k*R*T_in/(k - 1))*((Pc/Pe)^((k - 1)/k) - 1)   "The compression process is assumed to be isentropic";

      H_out = H_in + (k*R*T_in/(k - 1)) * ((Pc/Pe)^((k - 1)/k) - 1);

      connect(Power_Shaft, Power_Shaft)
        annotation (Line(points={{0,110},{0,110}}, color={0,0,127}));
     annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              lineThickness=1), Bitmap(
              extent={{-90,-88},{94,94}},
              imageSource="PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iaXNvLTg4NTktMSI/Pg0KPCEtLSBHZW5lcmF0b3I6IEFkb2JlIElsbHVzdHJhdG9yIDE5LjAuMCwgU1ZHIEV4cG9ydCBQbHVnLUluIC4gU1ZHIFZlcnNpb246IDYuMDAgQnVpbGQgMCkgIC0tPg0KPHN2ZyB2ZXJzaW9uPSIxLjEiIGlkPSJMYXllcl8xIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHhtbG5zOnhsaW5rPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hsaW5rIiB4PSIwcHgiIHk9IjBweCINCgkgdmlld0JveD0iMCAwIDUxMiA1MTIiIHN0eWxlPSJlbmFibGUtYmFja2dyb3VuZDpuZXcgMCAwIDUxMiA1MTI7IiB4bWw6c3BhY2U9InByZXNlcnZlIj4NCjxnPg0KCTxnPg0KCQk8cGF0aCBkPSJNNDkzLjczNiwyMDguMDc4Yy01LjU1LTEyLjEyLTE3LjYzNS0xOS45NTMtMzAuNzg4LTE5Ljk1M2gtNzAuNjY0di0wLjUzNGMwLTYuMzEtMy41NDctMTEuODAzLTguNzUtMTQuNg0KCQkJYzExLjA0OS0xMi4yNDQsMTcuMzAxLTI4LjE3NSwxNy4zMDEtNDUuMjU4YzAtMzcuNDI2LTMwLjQ0OC02Ny44NzUtNjcuODc1LTY3Ljg3NWMtOC45NTIsMC0xNi45NDgsMi4xNDItMTYuOTU0LDIuMTQzDQoJCQlsLTE1NC43ODIsNDEuNjQ5Yy0xOC44MzcsNC4yMzUtMzIuOTU3LDIxLjA4Ny0zMi45NTcsNDEuMTg0YzAsOS44ODcsMy40MjYsMTguOTg1LDkuMTQsMjYuMTg4aC0xLjEyMw0KCQkJYy05LjEzNiwwLTE2LjU2OCw3LjQzMi0xNi41NjgsMTYuNTY4djAuNTM0SDQ5LjA1MmMtMTMuMTU0LDAtMjUuMjM5LDcuODMzLTMwLjc4OCwxOS45NTRDNi4zMTUsMjM0LjE3NywwLDI2NS41MzIsMCwyOTguNzU2DQoJCQljMCwzMy4yMjMsNi4zMTYsNjQuNTc5LDE4LjI2NCw5MC42NzdjNS41NSwxMi4xMiwxNy42MzUsMTkuOTUzLDMwLjc4OCwxOS45NTNoNDkuODkybC0zLjgxNywyNi43MjJoLTEuNTk5DQoJCQljLTQuNDI3LDAtOC4wMTcsMy41ODktOC4wMTcsOC4wMTdjMCw0LjQyNywzLjU4OSw4LjAxNyw4LjAxNyw4LjAxN2g1OS44NThjNC40MjcsMCw4LjAxNy0zLjU4OSw4LjAxNy04LjAxNw0KCQkJYzAtNC40MjctMy41ODktOC4wMTctOC4wMTctOC4wMTdoLTEuNTk4bC0zLjgxOC0yNi43MjJoMjE2LjA1OWwtMy44MTcsMjYuNzIyaC0xLjU5OGMtNC40MjcsMC04LjAxNywzLjU4OS04LjAxNyw4LjAxNw0KCQkJYzAsNC40MjcsMy41ODksOC4wMTcsOC4wMTcsOC4wMTdoNTkuODU4YzQuNDI3LDAsOC4wMTctMy41ODksOC4wMTctOC4wMTdjMC00LjQyNy0zLjU4OS04LjAxNy04LjAxNy04LjAxN2gtMS41OThsLTMuODE4LTI2LjcyMg0KCQkJaDQ5Ljg5M2MxMy4xNTQsMCwyNS4yMzktNy44MzMsMzAuNzg4LTE5Ljk1NEM1MDUuNjg1LDM2My4zMzQsNTEyLDMzMS45NzksNTEyLDI5OC43NTYNCgkJCUM1MTIsMjY1LjUzMiw1MDUuNjg0LDIzNC4xNzYsNDkzLjczNiwyMDguMDc4eiBNMzMyLjk2LDc1Ljg5MWMyOC41ODUsMCw1MS44NDEsMjMuMjU2LDUxLjg0MSw1MS44NDENCgkJCWMwLDE1LjM5Ni02LjY3NSwyOS41NTktMTguMTcxLDM5LjMyNnYtNDcuODc3YzAtOS4xMzYtNy40MzItMTYuNTY4LTE2LjU2OC0xNi41NjhoLTM0LjIwNWMtOS4xMzYsMC0xNi41NjgsNy40MzItMTYuNTY4LDE2LjU2OA0KCQkJdjQ3Ljg3N2MtMTEuNDk2LTkuNzY4LTE4LjE3MS0yMy45MjktMTguMTcxLTM5LjMyNkMyODEuMTE5LDk5LjE0NywzMDQuMzc1LDc1Ljg5MSwzMzIuOTYsNzUuODkxeiBNMzE1LjMyNCwxMzYuODE4di0xNy42MzcNCgkJCWMwLTAuMjk1LDAuMjM5LTAuNTM0LDAuNTM0LTAuNTM0aDM0LjIwNWMwLjI5NSwwLDAuNTM0LDAuMjM5LDAuNTM0LDAuNTM0djE3LjYzN0gzMTUuMzI0eiBNMzUwLjU5NywxNTIuODUydjE4LjE3MWgtMzUuMjczDQoJCQl2LTE4LjE3MUgzNTAuNTk3eiBNMjc3LjI0Niw4OS4wMzJjLTcuNjU1LDEwLjk4Ny0xMi4xNjEsMjQuMzI0LTEyLjE2MSwzOC43YzAsMTYuMTgyLDUuNjEzLDMxLjMyOSwxNS41OTYsNDMuMjloLTc3LjExMQ0KCQkJYzUuNzE0LTcuMjAzLDkuMTQtMTYuMzAxLDkuMTQtMjYuMTg4Yy0wLjAwMS0xMy44NDUtNi43MDEtMjYuMTU1LTE3LjAyOS0zMy44NTlMMjc3LjI0Niw4OS4wMzJ6IE0xNzAuNDg5LDExOC42NDcNCgkJCWMxNC40NCwwLDI2LjE4OCwxMS43NDgsMjYuMTg4LDI2LjE4OHMtMTEuNzQ4LDI2LjE4OC0yNi4xODgsMjYuMTg4cy0yNi4xODgtMTEuNzQ4LTI2LjE4OC0yNi4xODgNCgkJCVMxNTYuMDQ5LDExOC42NDcsMTcwLjQ4OSwxMTguNjQ3eiBNMTM1Ljc0OSwxODcuNTkxYzAtMC4yOTUsMC4yMzktMC41MzQsMC41MzQtMC41MzRoMjM5LjQzMmMwLjI5NSwwLDAuNTM0LDAuMjM5LDAuNTM0LDAuNTM0DQoJCQl2MTcuMTAyYzAsMC4yOTUtMC4yMzksMC41MzQtMC41MzQsMC41MzRIMTM2LjI4NGMtMC4yOTUsMC0wLjUzNC0wLjIzOS0wLjUzNC0wLjUzNFYxODcuNTkxeiBNMTExLjMyMyw0MzYuMTA5bDYuMTk1LTQzLjM2Ng0KCQkJYzAuMDM3LTAuMjYyLDAuMjY1LTAuNDU4LDAuNTI5LTAuNDU4aDEwLjgyMWMwLjI2NCwwLDAuNDkyLDAuMTk4LDAuNTI5LDAuNDU4bDYuMTkzLDQzLjM2NkgxMTEuMzIzeiBNMzc2LjQwOSw0MzYuMTA5DQoJCQlsNi4xOTUtNDMuMzY2YzAuMDM3LTAuMjYyLDAuMjY1LTAuNDU4LDAuNTI5LTAuNDU4aDEwLjgyMWMwLjI2NCwwLDAuNDkyLDAuMTk4LDAuNTI5LDAuNDU4bDYuMTkzLDQzLjM2NkgzNzYuNDA5eg0KCQkJIE00NzkuMTU4LDM4Mi43NThjLTIuOTQ2LDYuNDM2LTkuMzA5LDEwLjU5NS0xNi4yMSwxMC41OTVoLTIuMjU1VjI2NC41NTFjMC00LjQyNy0zLjU4OS04LjAxNy04LjAxNy04LjAxNw0KCQkJYy00LjQyNywwLTguMDE3LDMuNTg5LTguMDE3LDguMDE3djEyOC44MDJoLTMzLjg5NWwtMC40MTItMi44NzdjLTEuMTU5LTguMTEtOC4yMDktMTQuMjI1LTE2LjQwMS0xNC4yMjVoLTEwLjgyMQ0KCQkJYy04LjE5MiwwLTE1LjI0Miw2LjExNS0xNi40MDEsMTQuMjI1bC0wLjQxLDIuODc3SDE0NS42OGwtMC40MTEtMi44NzdjLTEuMTU5LTguMTEtOC4yMDktMTQuMjI1LTE2LjQwMS0xNC4yMjVoLTEwLjgyMQ0KCQkJYy04LjE5MiwwLTE1LjI0Miw2LjExNS0xNi40MDEsMTQuMjI1bC0wLjQxLDIuODc3SDY3LjM0VjI2NC41NTFjMC00LjQyNy0zLjU4OS04LjAxNy04LjAxNy04LjAxN3MtOC4wMTcsMy41ODktOC4wMTcsOC4wMTcNCgkJCXYxMjguODAyaC0yLjI1NGMtNi45MDEsMC0xMy4yNjQtNC4xNTktMTYuMjEtMTAuNTk1Yy0xMC44NC0yMy42NzUtMTYuODA5LTUzLjUwOC0xNi44MDktODQuMDAyczUuOTY5LTYwLjMyNiwxNi44MDgtODQuMDAyDQoJCQljMi45NDYtNi40MzYsOS4zMDktMTAuNTk1LDE2LjIxLTEwLjU5NWgyLjI1NXYyNi4xODhjMCw0LjQyNywzLjU4OSw4LjAxNyw4LjAxNyw4LjAxN3M4LjAxNy0zLjU4OSw4LjAxNy04LjAxN3YtMjYuMTg4aDUyLjM3Ng0KCQkJdjAuNTM0YzAsOS4xMzYsNy40MzIsMTYuNTY4LDE2LjU2OCwxNi41NjhoMTcuNjM3djAuNTM0YzAsNC40MjcsMy41ODksOC4wMTcsOC4wMTcsOC4wMTdzOC4wMTctMy41ODksOC4wMTctOC4wMTd2LTAuNTM0aDE3Mi4wOTINCgkJCXYwLjUzNGMwLDQuNDI3LDMuNTg5LDguMDE3LDguMDE3LDguMDE3czguMDE3LTMuNTg5LDguMDE3LTguMDE3di0wLjUzNGgxNy42MzdjOS4xMzYsMCwxNi41NjgtNy40MzIsMTYuNTY4LTE2LjU2OHYtMC41MzRoNTIuMzc2DQoJCQl2MjYuMTg4YzAsNC40MjcsMy41ODksOC4wMTcsOC4wMTcsOC4wMTdjNC40MjcsMCw4LjAxNy0zLjU4OSw4LjAxNy04LjAxN3YtMjYuMTg4aDIuMjU0YzYuOTAxLDAsMTMuMjY0LDQuMTU5LDE2LjIxLDEwLjU5NQ0KCQkJYzEwLjg0LDIzLjY3NSwxNi44MDksNTMuNTA4LDE2LjgwOSw4NC4wMDJTNDg5Ljk5OCwzNTkuMDgyLDQ3OS4xNTgsMzgyLjc1OHoiLz4NCgk8L2c+DQo8L2c+DQo8Zz4NCjwvZz4NCjxnPg0KPC9nPg0KPGc+DQo8L2c+DQo8Zz4NCjwvZz4NCjxnPg0KPC9nPg0KPGc+DQo8L2c+DQo8Zz4NCjwvZz4NCjxnPg0KPC9nPg0KPGc+DQo8L2c+DQo8Zz4NCjwvZz4NCjxnPg0KPC9nPg0KPGc+DQo8L2c+DQo8Zz4NCjwvZz4NCjxnPg0KPC9nPg0KPGc+DQo8L2c+DQo8L3N2Zz4NCg==",
              fileName="//nas.ads.mwn.de/ge43kok/Desktop/SVGs/air-compressor-svgrepo-com.svg")}),
                                                                    Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StopTime=1000,
          __Dymola_NumberOfIntervals=10000,
          __Dymola_Algorithm="Dassl"));
    end Compressor;

  end Heat_Pump;

  package Experiment

    model Example
      Heat_Pump.Compressor compressor(V_c=50e-5)  annotation (Placement(
            transformation(
            extent={{10,10},{-10,-10}},
            rotation=270,
            origin={-56,12})));
      Heat_Pump.ExpanValve expanValve        annotation (Placement(
            transformation(
            extent={{10,10},{-10,-10}},
            rotation=90,
            origin={72,2})));
      Modelica.Blocks.Sources.RealExpression WaterTemperature_K(y=293.15)
        annotation (Placement(transformation(extent={{-48,42},{-34,62}})));
      Modelica.Blocks.Sources.RealExpression WaterMassFlow_kg_per_s(y=0.15)
        annotation (Placement(transformation(extent={{-36,30},{-22,50}})));
      Modelica.Blocks.Sources.RealExpression AirTemperature_K(y=283.15)
        annotation (Placement(transformation(extent={{68,-62},{54,-42}})));
      Heat_Pump.Condensor_R134a condensor_R134a(
        n=1,
        Ro=12.5e-3,
        Ri=11.5e-3,
        Alpa_i_vapor=400,
        Alpa_i_satu=2000,
        Alpa_o=2000,
        Cp_wall=385,
        D_wall=8960,
        Volume_water=0.005,
        Pc_initial=600000,
        H_refri_out_initial=320e3)
        annotation (Placement(transformation(extent={{-6,54},{14,74}})));
      Heat_Pump.Evaporator_R134a evaporator_R134a(
        n=1,
        Ro=9.5e-3,
        Ri=8.25e-3,
        Alpa_i_vapor=400,
        Alpa_i_satu=2000,
        Alpa_o=300,
        Cp_wall=385,
        D_wall=8960,
        Pe_initial=400000,
        H_refri_out_initial=410e3,
        L1_initial=5)
        annotation (Placement(transformation(extent={{16,-72},{-4,-52}})));
      Modelica.Blocks.Sources.RealExpression Speed_rads_per_s(y=60)
        annotation (Placement(transformation(extent={{-24,2},{-38,22}})));
      Modelica.Blocks.Math.Feedback feedback
        annotation (Placement(transformation(extent={{-2,-36},{14,-20}})));
      Modelica.Blocks.Continuous.PI PI(k=-1e-6, T=5)
        annotation (Placement(transformation(extent={{24,2},{36,14}})));
      Modelica.Blocks.Math.Add add
        annotation (Placement(transformation(extent={{44,-4},{56,8}})));
      Modelica.Blocks.Sources.Constant BaseOrificeArea(k=0.5e-6)
        annotation (Placement(transformation(extent={{20,-52},{32,-40}})));
      Modelica.Blocks.Sources.RealExpression SuperHeatRef(y=15)
        annotation (Placement(transformation(extent={{8,-22},{-6,-2}})));
    equation
      connect(compressor.Pc, expanValve.Pc) annotation (Line(points={{-48,1},{
              -48,-2},{-18,-2},{-18,32},{70,32},{70,13}}, color={0,0,127}));
      connect(evaporator_R134a.Pe, expanValve.Pe) annotation (Line(points={{-5,
              -54},{-16,-54},{-16,20},{66,20},{66,13}}, color={0,0,127}));
      connect(condensor_R134a.D_refri_out, expanValve.D_in)
        annotation (Line(points={{15,66},{74,66},{74,13}}, color={0,0,127}));
      connect(condensor_R134a.H_refri_out, expanValve.H_in)
        annotation (Line(points={{15,70},{78,70},{78,13}}, color={0,0,127}));
      connect(evaporator_R134a.T_air, AirTemperature_K.y) annotation (Line(
            points={{17,-56},{50,-56},{50,-52},{53.3,-52}}, color={0,0,127}));
      connect(evaporator_R134a.H_refri_in, expanValve.H_out)
        annotation (Line(points={{17,-60},{74,-60},{74,-9}}, color={0,0,127}));
      connect(evaporator_R134a.MassFlow_refri_in, expanValve.MassFlow_out)
        annotation (Line(points={{17,-64},{78,-64},{78,-9}}, color={0,0,127}));
      connect(compressor.MassFlow, evaporator_R134a.MassFlow_refri_out)
        annotation (Line(points={{-62,23},{-62,28},{-74,28},{-74,-78},{22,-78},
              {22,-68},{17,-68}}, color={0,0,127}));
      connect(evaporator_R134a.H_refri_out, compressor.H_in) annotation (Line(
            points={{-5,-58},{-56,-58},{-56,1}}, color={0,0,127}));
      connect(evaporator_R134a.D_refri_out, compressor.D_in) annotation (Line(
            points={{-5,-62},{-60,-62},{-60,1}}, color={0,0,127}));
      connect(evaporator_R134a.T_refri_out, compressor.T_in) annotation (Line(
            points={{-5,-66},{-64,-66},{-64,1}}, color={0,0,127}));
      connect(compressor.Pe, expanValve.Pe) annotation (Line(points={{-52,1},{
              -52,-54},{-16,-54},{-16,20},{66,20},{66,13}}, color={0,0,127}));
      connect(compressor.H_out, condensor_R134a.H_refri_in)
        annotation (Line(points={{-58,23},{-58,64},{-7,64}}, color={0,0,127}));
      connect(condensor_R134a.MassFlow_refri_in, evaporator_R134a.MassFlow_refri_out)
        annotation (Line(points={{-7,68},{-62,68},{-62,28},{-74,28},{-74,-78},{
              22,-78},{22,-68},{17,-68}}, color={0,0,127}));
      connect(condensor_R134a.T_water_in, WaterTemperature_K.y) annotation (
          Line(points={{-7,60},{-26,60},{-26,52},{-33.3,52}}, color={0,0,127}));
      connect(condensor_R134a.MassFlow_water_in, WaterMassFlow_kg_per_s.y)
        annotation (Line(points={{-7,56},{-16,56},{-16,40},{-21.3,40}}, color={
              0,0,127}));
      connect(condensor_R134a.MassFlow_refri_out, expanValve.MassFlow_out)
        annotation (Line(points={{-7,72},{-24,72},{-24,80},{88,80},{88,-18},{78,
              -18},{78,-9}}, color={0,0,127}));
      connect(condensor_R134a.Pc, expanValve.Pc)
        annotation (Line(points={{15,62},{70,62},{70,13}}, color={0,0,127}));
      connect(Speed_rads_per_s.y, compressor.Speed)
        annotation (Line(points={{-38.7,12},{-45,12}}, color={0,0,127}));
      connect(feedback.u2, evaporator_R134a.Superheat)
        annotation (Line(points={{6,-34.4},{6,-51}}, color={0,0,127}));
      connect(add.u1, PI.y) annotation (Line(points={{42.8,5.6},{42,5.6},{42,6},
              {40,6},{40,8},{36.6,8}}, color={0,0,127}));
      connect(add.y, expanValve.A_v)
        annotation (Line(points={{56.6,2},{61,2}}, color={0,0,127}));
      connect(BaseOrificeArea.y, add.u2) annotation (Line(points={{32.6,-46},{
              40,-46},{40,-1.6},{42.8,-1.6}}, color={0,0,127}));
      connect(feedback.y, PI.u) annotation (Line(points={{13.2,-28},{16,-28},{
              16,8},{22.8,8}}, color={0,0,127}));
      connect(SuperHeatRef.y, feedback.u1) annotation (Line(points={{-6.7,-12},
              {-8,-12},{-8,-28},{-0.4,-28}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Example;
  end Experiment;

  annotation (
    uses(Modelica(version = "3.2.3"),
      ElectrifiedPowertrains(version="1.3.1"),
      DassaultSystemes(version="1.3.1"),
      FluidPower(version="2019.2"),
      ClaytexFluid(version="2019.2"),
      Cooling(version="1.3.2"),
      HVACLib(version="2.8.0"),
      ElectricPowerSystems(version="1.3"),
      DymolaModels(version="1.0.1"),
      Battery(version="2.1.3"),
      ThermalSystems(version="1.5.0"),
      TSMedia(version="1.5.0"),
      Modelica_LinearSystems2(version="2.3.5")), Icon(graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          lineThickness=1),
        Line(
          points={{100,100},{-100,-100}},
          color={0,0,0},
          thickness=1),
        Text(
          extent={{-80,72},{-30,6}},
          lineColor={244,125,35},
          lineThickness=1,
          textString="E"),
        Text(
          extent={{28,-8},{78,-74}},
          lineColor={238,46,47},
          lineThickness=1,
          textString="H")}));

end MB_Library_HeatPump;
