namespace tRNA;

internal sealed class Program
{
	static void Main(string[] args)
	{
		for (int i = 0; i < args.Length; i++)
		{
			Console.WriteLine($"Argument {i}: {args[i]}");
		}
	}
}
