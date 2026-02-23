import {
  Bar,
  BarChart,
  CartesianGrid,
  Legend,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis
} from "recharts";

export default function PortfolioGreekBarChart({ chartData, rows, colors }) {
  return (
    <ResponsiveContainer width="100%" height="100%">
      <BarChart data={chartData}>
        <CartesianGrid strokeDasharray="3 3" />
        <XAxis dataKey="greek" />
        <YAxis />
        <Tooltip formatter={(value) => Number(value).toFixed(3)} />
        <Legend />
        {rows.map((row, index) => (
          <Bar key={row.id} dataKey={`leg${index}`} stackId="a" fill={colors[index % colors.length]} name={row.label} />
        ))}
      </BarChart>
    </ResponsiveContainer>
  );
}
