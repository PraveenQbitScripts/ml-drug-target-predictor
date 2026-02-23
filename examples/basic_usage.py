#!/usr/bin/env python3
"""
Basic usage example for ML Drug Target Predictor
"""

import sys
import os
sys.path.append('..')

def main():
    """Example usage of ML Drug Target Predictor"""
    
    # Example protein sequences (rice genes)
    proteins = [
        "Os01g0100100",
        "Os01g0100200", 
        "Os01g0100300"
    ]
    
    print("🤖 ML Drug Target Predictor - Basic Usage Example")
    print("=" * 50)
    print(f"Analyzing {len(proteins)} proteins...")
    
    # This would normally call the actual predictor
    # from Script_10_Unified_Improved_1 import DrugTargetPredictor
    # predictor = DrugTargetPredictor()
    # results = predictor.predict_targets(proteins)
    
    # Simulated results for demonstration
    simulated_results = {
        'total_proteins': len(proteins),
        'high_confidence_targets': 2,
        'medium_confidence_targets': 1,
        'processing_time': '45 seconds'
    }
    
    print(f"\n📊 Results:")
    print(f"  Total proteins analyzed: {simulated_results['total_proteins']}")
    print(f"  High-confidence targets: {simulated_results['high_confidence_targets']}")
    print(f"  Medium-confidence targets: {simulated_results['medium_confidence_targets']}")
    print(f"  Processing time: {simulated_results['processing_time']}")
    
    print(f"\n✅ Analysis complete!")
    print("📄 Check output directory for detailed results")

if __name__ == "__main__":
    main()
