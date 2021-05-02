
class HeatTransferSurface
{
  public:
    HeatTransferSurface(const double area, const double r_value) :
      m_area(area),
      m_r_value(r_value)
  {}


    // given the delta_T (in Kelvin), computes the rate of heat transfer
    // (W, or J/s) through the surface
    double computeHeatRate(const double delta_T)
    {
      return delta_T * m_area / m_r_value;
    }

  private:
    double m_area;     // units: m^2
    double m_r_value;  // K m^2/W
};


class AirExchangeHeatTransfer
{
  public:
    AirExchangeHeatTransfer(const double volume, const double ach) :
      m_volume(volume),
      m_ach(ach)
  {}

    double computeHeatRate(const double delta_T)
    {
      // TODO: the 0.018 figure comes from
      // https://www.pexuniverse.com/calculate-heat-loss
      // should figure out the underlying physics
      // convert air changes per hour to air changes per second
      // This gives a result in Joules per second (or W)
      return m_volume * (m_ach / 3600) * m_rho * m_Cp * delta_T;
    }

  private:
    double m_volume;  // m^3
    double m_ach;     // number of air changes per hour
    double m_rho = 1.225;  // density of air, kg/m^3
    double m_Cp  = 1006;   // specific heat of air at constant pressure,
                           // J/(kg K)
};
